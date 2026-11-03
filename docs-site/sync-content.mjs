#!/usr/bin/env node
/**
 * Syncs ALL documentation content into Starlight's src/content/docs/ tree.
 *
 * ARCHITECTURAL RULE: src/content/docs/ is 100% GENERATED.
 * ─────────────────────────────────────────────────────────
 * Every page originates in docs/ or, for a study, in the README that sits next
 * to the study's code under benchmarks/.  Nothing is ever authored directly in
 * src/content/docs/.  The entire directory is .gitignored.
 *
 * Content categories:
 *   1. SITE pages  (docs/site/*.mdx)          — copied verbatim (already have frontmatter)
 *   2. SYNCED docs (docs/**\/*.md)             — H1 extracted as title, frontmatter added,
 *                                               image paths & links rewritten
 *   3. STUDIES     (benchmarks/**\/README.md)  — same transform, plus the study's
 *                                               figures/ copied to public/img/studies/<slug>/
 *                                               and its other relative links sent to GitHub
 *
 * Run automatically via `npm run build` / `npm run dev`.
 */

import { readFileSync, writeFileSync, mkdirSync, existsSync, cpSync, rmSync } from 'fs';
import { dirname, join, posix } from 'path';

const REPO = join(import.meta.dirname, '..');
const DOCS = join(REPO, 'docs');
const OUT  = join(import.meta.dirname, 'src', 'content', 'docs');
const PUB  = join(import.meta.dirname, 'public');
const BASE = '/mp-numerics';  // Astro base path
const GITHUB = 'https://github.com/stillwater-sc/mp-numerics';

// ── Site pages (MDX with Starlight components) ────────────────────
// These are copied verbatim — they already contain frontmatter.
const SITE_FILES = {
  'site/index.mdx': 'index.mdx',
};

// ── Synced docs (source path relative to docs/ → dest relative to content/docs/) ──
const FILE_MAP = {
  'getting-started/index.md': 'getting-started/index.md',
  'roadmap.md':               'roadmap.md',
};

// ── Studies (study directory relative to the repo root → page slug) ──
// The README in each directory is the page; its figures/ are copied next to
// the site's images, so the README renders the same on GitHub and here.
const STUDY_MAP = {
  'benchmarks/linalg/givens_householder': 'studies/linalg/givens-householder',
  'benchmarks/linalg/decompositions':     'studies/linalg/decompositions',
};

// ── Link lookup ───────────────────────────────────────────────────

function slugUrl(slug) {
  return `${BASE}/${slug.endsWith('/') ? slug : slug + '/'}`;
}

function buildLinkLookup() {
  const lookup = {};
  // Keys are paths relative to docs/, as a docs/ page resolves its links.
  for (const [src, dest] of Object.entries(FILE_MAP)) {
    lookup[src] = slugUrl(dest.replace(/\.md$/, '').replace(/(^|\/)index$/, '$1'));
  }
  for (const [src, dest] of Object.entries(SITE_FILES)) {
    lookup[src] = slugUrl(dest.replace(/\.mdx?$/, '').replace(/(^|\/)index$/, '$1'));
  }
  // A study README sits outside docs/, so a link to it normalizes with a
  // leading "../".
  for (const [dir, slug] of Object.entries(STUDY_MAP)) {
    lookup[`../${dir}/README.md`] = slugUrl(slug);
  }
  return lookup;
}

const LINK_LOOKUP = buildLinkLookup();

// ── Transforms ────────────────────────────────────────────────────

function rewriteLinks(content, srcRelative) {
  const srcDir = posix.dirname(srcRelative);
  // Match ](path.md) and ](path.md#anchor), or .mdx for a site page; the path
  // is resolved via the link lookup and a trailing #section anchor (if any) is
  // preserved verbatim, so cross-page section links survive the rewrite.
  return content.replace(/\]\(([^)#]+\.mdx?)(#[^)]*)?\)/g, (match, target, anchor = '') => {
    if (target.startsWith('http://') || target.startsWith('https://')) return match;
    const resolved = posix.normalize(posix.join(srcDir, target));
    const url = LINK_LOOKUP[resolved];
    return url ? `](${url}${anchor})` : match;
  });
}

function extractTitle(content) {
  const match = content.match(/^#\s+(.+)$/m);
  return match ? match[1].trim() : 'Untitled';
}

function stripFirstHeading(content) {
  return content.replace(/^#\s+.+\n*/m, '');
}

function rewriteImagePaths(content) {
  // docs/ pages use a relative "../img/..." path so the images also render in
  // GitHub's markdown view. Collapse any leading "../" to the site-absolute
  // "/mp-numerics/img/..." for the built site.
  return content.replace(/\]\((?:\.\.\/)*img\//g, `](${BASE}/img/`);
}

function rewriteStudyLinks(content, studyDir, slug) {
  // figures/ are served from the site; every other relative link that is not
  // a page (results/, sources, scripts) goes to the file on GitHub.
  return content.replace(/\]\(([^)#\s]+)(#[^)]*)?\)/g, (match, target, anchor = '') => {
    if (/^[a-z]+:/.test(target) || target.startsWith('/')) return match;
    if (/\.mdx?$/.test(target)) return match;
    const resolved = posix.normalize(posix.join(studyDir, target));
    const figures = `${studyDir}/figures/`;
    if (resolved.startsWith(figures)) {
      return `](${BASE}/img/studies/${slug.split('/').pop()}/${resolved.slice(figures.length)}${anchor})`;
    }
    return `](${GITHUB}/blob/main/${resolved}${anchor})`;
  });
}

function frontmatter(title, extra = '') {
  const safeTitle = title.replace(/\\/g, '\\\\').replace(/"/g, '\\"');
  return `---\ntitle: "${safeTitle}"\n${extra}---\n\n`;
}

function addFrontmatter(content, srcRelative) {
  let body = stripFirstHeading(content);
  body = rewriteImagePaths(body);
  body = rewriteLinks(body, srcRelative);
  return frontmatter(extractTitle(content)) + body;
}

function studyPage(content, studyDir, slug) {
  let body = stripFirstHeading(content);
  body = rewriteStudyLinks(body, studyDir, slug);
  body = rewriteLinks(body, `../${studyDir}/README.md`);
  // The edit link points at the README, not at docs/.
  const edit = `editUrl: "${GITHUB}/edit/main/${studyDir}/README.md"\n`;
  const source =
    `:::note[Source]\nThis page is \`${studyDir}/README.md\`. The code, the raw ` +
    `results (CSV with \`.sysinfo\` sidecars) and the plotting scripts are in ` +
    `[\`${studyDir}/\`](${GITHUB}/tree/main/${studyDir}).\n:::\n\n`;
  return frontmatter(extractTitle(content), edit) + source + body;
}

// ── File operations ───────────────────────────────────────────────

function writeOut(destRelative, content) {
  const destPath = join(OUT, destRelative);
  mkdirSync(dirname(destPath), { recursive: true });
  writeFileSync(destPath, content);
}

function readSource(srcPath) {
  if (!existsSync(srcPath)) {
    console.warn(`  SKIP (not found): ${srcPath}`);
    return null;
  }
  return readFileSync(srcPath, 'utf-8');
}

// ── Main ──────────────────────────────────────────────────────────

// Clear stale Astro data store cache
const astroCache = join(import.meta.dirname, 'node_modules', '.astro');
if (existsSync(astroCache)) {
  rmSync(astroCache, { recursive: true });
}

// Wipe the generated directories — they are 100% generated
for (const dir of [OUT, join(PUB, 'img')]) {
  if (existsSync(dir)) rmSync(dir, { recursive: true });
}
mkdirSync(OUT, { recursive: true });

console.log('Syncing docs/ and benchmarks/ → docs-site/src/content/docs/ ...');

// 1. Copy site pages (MDX with components, already have frontmatter)
for (const [src, dest] of Object.entries(SITE_FILES)) {
  const content = readSource(join(DOCS, src));
  if (content === null) continue;
  writeOut(dest, content);
  console.log(`  site: ${src} → ${dest}`);
}

// 2. Sync docs/ markdown files (add frontmatter, rewrite links)
for (const [src, dest] of Object.entries(FILE_MAP)) {
  const content = readSource(join(DOCS, src));
  if (content === null) continue;
  writeOut(dest, addFrontmatter(content, src));
  console.log(`  sync: ${src} → ${dest}`);
}

// 3. Sync the studies and copy their figures
for (const [dir, slug] of Object.entries(STUDY_MAP)) {
  const content = readSource(join(REPO, dir, 'README.md'));
  if (content === null) continue;
  writeOut(`${slug}.md`, studyPage(content, dir, slug));
  console.log(`  study: ${dir}/README.md → ${slug}.md`);
  const figures = join(REPO, dir, 'figures');
  if (existsSync(figures)) {
    const dest = join(PUB, 'img', 'studies', slug.split('/').pop());
    mkdirSync(dest, { recursive: true });
    cpSync(figures, dest, { recursive: true });
    console.log(`         ${dir}/figures/ → public/img/studies/${slug.split('/').pop()}/`);
  }
}

// 4. Copy docs/img/ to public/ (served as static assets at /mp-numerics/img/)
const imgSrc = join(DOCS, 'img');
if (existsSync(imgSrc)) {
  const imgDest = join(PUB, 'img');
  mkdirSync(imgDest, { recursive: true });
  cpSync(imgSrc, imgDest, { recursive: true });
  console.log('  Copied docs/img/ → public/img/');
}

console.log('Done.');

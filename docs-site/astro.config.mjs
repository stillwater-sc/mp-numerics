import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';

export default defineConfig({
  site: 'https://stillwater-sc.github.io',
  base: '/mp-numerics',
  integrations: [
    starlight({
      title: 'mp-numerics -- Mixed-Precision Numerical Methods',
      description:
        'Classic numerical methods studied across number systems, composing MTL5 and Universal',
      social: [
        {
          icon: 'github',
          label: 'GitHub',
          href: 'https://github.com/stillwater-sc/mp-numerics',
        },
      ],
      editLink: {
        // Turns the edit link on. Every page sets its own editUrl (from
        // sync-content.mjs, or in a docs/site/*.mdx frontmatter), because
        // this base would be joined with the generated src/content/docs/ path.
        baseUrl:
          'https://github.com/stillwater-sc/mp-numerics/edit/main/docs/',
      },
      customCss: ['./src/styles/custom.css'],
      sidebar: [
        {
          label: 'Getting Started',
          autogenerate: { directory: 'getting-started' },
        },
        {
          label: 'Studies',
          items: [
            {
              label: 'Linear algebra',
              items: [
                { slug: 'studies/linalg/givens-householder' },
                { slug: 'studies/linalg/decompositions' },
              ],
            },
          ],
        },
        { slug: 'roadmap' },
        {
          label: 'Related',
          items: [
            {
              label: 'MTL5 documentation',
              link: 'https://stillwater-sc.github.io/mtl5/',
            },
            {
              label: 'Givens/Householder explorer (MTL5)',
              link: 'https://stillwater-sc.github.io/mtl5/algorithms/givens-householder-explorer/',
            },
            {
              label: 'Universal documentation',
              link: 'https://stillwater-sc.github.io/universal/',
            },
          ],
        },
      ],
    }),
  ],
});

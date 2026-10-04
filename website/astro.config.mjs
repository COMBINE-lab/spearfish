// @ts-check
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';
import remarkGfm from 'remark-gfm';

// GitHub Pages: site + base. The site is meant for
// https://combine-lab.github.io/spearfish/ (not deployed yet).
export default defineConfig({
  site: 'https://combine-lab.github.io',
  base: '/spearfish',
  // GitHub-Flavored Markdown is applied to `.md` by default but not to `.mdx`;
  // add remark-gfm explicitly so tables in the `.mdx` pages render.
  markdown: {
    remarkPlugins: [remarkGfm],
  },
  integrations: [
    starlight({
      title: 'spearfish',
      description:
        'spearfish: fast long-read transcriptome mapping by delta mapping. Every compatible isoform, scored, with each piece of evidence computed once.',
      logo: {
        light: './src/assets/spearfish-logo.svg',
        dark: './src/assets/spearfish-logo-dark.svg',
        alt: 'spearfish',
      },
      social: [
        {
          icon: 'github',
          label: 'GitHub',
          href: 'https://github.com/COMBINE-lab/spearfish',
        },
      ],
      editLink: {
        baseUrl: 'https://github.com/COMBINE-lab/spearfish/edit/main/website/',
      },
      sidebar: [
        {
          label: 'Getting started',
          items: [
            { label: 'Introduction', slug: 'getting-started/introduction' },
            { label: 'Installation', slug: 'getting-started/installation' },
            { label: 'Quick start', slug: 'getting-started/quick-start' },
          ],
        },
        {
          label: 'Guide',
          items: [
            { label: 'Building an index', slug: 'guide/indexing' },
            { label: 'Mapping reads', slug: 'guide/mapping' },
            { label: 'Output format', slug: 'guide/output' },
            { label: 'Quantification', slug: 'guide/quantification' },
          ],
        },
        {
          label: 'How it works',
          items: [
            { label: 'Delta mapping', slug: 'concepts/delta-mapping' },
            { label: 'Inside spearfish', slug: 'concepts/inside' },
          ],
        },
        {
          label: 'About',
          items: [
            { label: 'Citation', slug: 'about/citation' },
            { label: 'License', slug: 'about/license' },
          ],
        },
      ],
    }),
  ],
});

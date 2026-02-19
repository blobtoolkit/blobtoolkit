# BlobToolKit Documentation Site

This is a Jekyll-based documentation site for BlobToolKit. The site is built with Jekyll 4.4.1 and uses a custom theme based on jekyll-theme-minimal.

## Prerequisites

- Ruby 3.2.0 or later
- Bundler
- rbenv (recommended for Ruby version management)

## Setup

### 1. Install Ruby (if needed)

If you don't have Ruby 3.2.0 installed, use rbenv:

```bash
# Install rbenv if you haven't already
brew install rbenv

# Install Ruby 3.2.0
rbenv install 3.2.0

# Set it as the local version for this project
rbenv local 3.2.0
```

### 2. Install Dependencies

From the `docs` directory, install the required gems:

```bash
bundle install
```

This installs Jekyll and other dependencies specified in the `Gemfile`.

## Running Locally

### Quick Start (Recommended)

Use the development helper script that automatically generates navigation:

```bash
bash scripts/dev.sh
```

This script:

1. Generates navigation from folder structure
2. Installs dependencies if needed
3. Starts the Jekyll server with incremental builds

The site will be available at `http://localhost:4000`

### Manual Start

If you prefer to manage the process manually:

1. Generate navigation:

   ```bash
   ruby scripts/generate_nav.rb
   ```

2. Start Jekyll:
   ```bash
   bundle exec jekyll serve
   ```

### Options

- `--incremental` - Only regenerate changed files (faster)
- `--host localhost --port 4000` - Specify custom host/port
- `--draft` - Include draft posts

Example:

```bash
bundle exec jekyll serve --incremental
```

## Project Structure

```
docs/
├── index.md                    # Home page
├── installation/
│   └── index.md               # Installation overview
├── pipeline/
│   ├── index.md               # Pipeline overview
│   ├── getting-started.md
│   └── configuration.md
├── command-line/
│   ├── index.md               # Command-line overview
│   ├── commands.md
│   └── tutorials.md
├── viewer/
│   ├── index.md               # Viewer overview
│   ├── local-hosting.md
│   └── visualisations.md
├── specification/
│   ├── index.md               # Specification overview
│   ├── blobdir-format.md
│   └── validator.md
├── _layouts/
│   └── default.html           # Master layout template
├── assets/
│   ├── css/
│   │   └── style.scss         # Custom SCSS with branding
│   └── img/
│       └── logo.png           # Site logo
├── _config.yml                # Jekyll configuration
└── Gemfile                    # Ruby dependencies
```

## Making Changes

## Making Changes

### Adding/Editing Pages

1. Create or edit `.md` files in the appropriate directories
2. Add Jekyll front matter at the top of each file:

   ```yaml
   ---
   layout: default
   title: Page Title
   ---
   Page content here...
   ```

3. **Update navigation** by running the navigation generator script:

   ```bash
   ruby scripts/generate_nav.rb
   ```

   This script automatically:
   - Scans all folders in the docs directory
   - Finds all markdown files in each folder
   - Extracts titles from front matter
   - Generates `_data/navigation.yml`
   - The sidebar automatically updates based on the folder structure

4. Jekyll will regenerate the site when you save

### Understanding Dynamic Navigation

The navigation system uses a two-step approach:

1. **Build time**: `scripts/generate_nav.rb` scans the docs directory and generates `_data/navigation.yml`
2. **Render time**: The layout uses `site.data.navigation` to dynamically display:
   - On the home page: All top-level sections
   - On section pages: All pages within that section with the current page highlighted

**Important**: Run `ruby scripts/generate_nav.rb` whenever you:

- Add a new folder (section)
- Add a new markdown file to a section
- Change a page title in the front matter

**Automatic updates**: You can run this script as part of your development workflow or add it to a pre-commit hook.

### Adding New Sections

To create a new documentation section:

1. Create a new folder in `/docs/`:

   ```bash
   mkdir docs/my-new-section
   ```

2. Create `docs/my-new-section/index.md` with front matter:

   ```yaml
   ---
   layout: default
   title: My New Section
   ---
   ```

3. Add other `.md` files to this folder as needed

4. Run the navigation generator:

   ```bash
   ruby scripts/generate_nav.rb
   ```

5. The new section will automatically appear in:
   - The sidebar on the home page
   - The main navigation menu (if you update `_config.yml`)

### Updating Styles

Edit `assets/css/style.scss`:

- SCSS variables and custom styles are at the top
- Media queries for responsive design use a 1100px breakpoint
- The Jekyll theme base styles are imported at the beginning

### Custom Theme

The site uses custom branding with:

- **Primary accent color**: `#9c528b` (purple)
- **Dark backgrounds**: `#00102E` and `#233452`
- **Fonts**: Comfortaa (headings), Open Sans (body)
- **Full-width header/footer** with fixed positioning
- **Responsive sidebar navigation** that collapses on mobile

## Responsive Design

- **Desktop (>1100px)**: Sidebar always visible on left
- **Mobile (<1100px)**: Sidebar hidden by default, toggle with hamburger menu button

The responsive behavior is controlled by CSS media queries in `style.scss`.

## Deployment

The site is deployed to GitHub Pages when changes are pushed to the repository. The `_site/` folder contains the compiled static HTML and is automatically generated by Jekyll.

## Troubleshooting

### Jekyll won't start

- Ensure you're using Ruby 3.2.0: `rbenv local 3.2.0`
- Run `bundle install` to ensure all gems are installed
- Check for syntax errors in `_config.yml`

### Changes not appearing

- Clear Jekyll cache: `rm -rf _site .jekyll-cache`
- Stop and restart the server with `bundle exec jekyll serve`

### CSS not updating

- Ensure `assets/css/style.scss` has proper SCSS syntax
- The file must start with `---` and `---` (Jekyll front matter) for SCSS compilation
- Restart the server after CSS changes

### Port already in use

- Use a different port: `bundle exec jekyll serve --port 4001`
- Or kill the process: `lsof -ti:4000 | xargs kill -9`

## Contributing

When adding new documentation:

1. Follow the existing directory structure
2. Use markdown formatting
3. Include front matter with `layout: default` and a `title`
4. Test locally before committing
5. Ensure links to internal pages don't include `.md` extension (Jekyll auto-converts)

## Resources

- [Jekyll Documentation](https://jekyllrb.com/docs/)
- [Markdown Guide](https://www.markdownguide.org/)
- [jekyll-theme-minimal](https://github.com/pages-themes/minimal)

#!/bin/bash
# Development helper script for the BlobToolKit docs site
# Generates navigation and starts Jekyll dev server

set -e

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOCS_DIR="$(dirname "$SCRIPT_DIR")"

echo "🔧 BlobToolKit Docs Development Setup"
echo "======================================="
echo ""

# Check if we're in the right directory
if [ ! -f "$DOCS_DIR/Gemfile" ]; then
    echo "❌ Error: Gemfile not found at $DOCS_DIR/Gemfile"
    exit 1
fi

cd "$DOCS_DIR"

# Clear Jekyll cache for fresh build
echo "🧹 Clearing Jekyll cache..."
rm -rf _site .jekyll-cache .jekyll-metadata
echo ""

# Generate navigation
echo "📋 Generating navigation from folder structure..."
ruby "$SCRIPT_DIR/generate_nav.rb"
echo ""

# Optional: offer to bundle install if gems aren't installed
if ! bundle check > /dev/null 2>&1; then
    echo "📦 Installing dependencies..."
    bundle install
    echo ""
fi

# Start Jekyll
echo "🚀 Starting Jekyll development server..."
echo "📍 Site will be available at http://localhost:4000"
echo ""
echo "💡 Tips:"
echo "  - Add new markdown files to any section folder"
echo "  - Update page titles in the front matter"
echo "  - Run 'ruby scripts/generate_nav.rb' to update navigation"
echo "  - Press Ctrl+C to stop the server"
echo ""

bundle exec jekyll serve --incremental

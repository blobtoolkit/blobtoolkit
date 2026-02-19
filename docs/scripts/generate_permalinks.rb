#!/usr/bin/env ruby

require 'fileutils'

# Automatically add permalinks to markdown files that need them
class PermalinkGenerator
  def initialize(docs_dir)
    @docs_dir = docs_dir
  end

  def generate
    Dir.glob(File.join(@docs_dir, '**/*.md')).each do |file|
      # Skip index.md files and files in hidden/special directories
      next if File.basename(file) == 'index.md'
      next if file.include?('/_')
      next if file.include?('/.git')
      
      process_file(file)
    end
  end

  private

  def process_file(file)
    content = File.read(file)
    
    # Parse frontmatter
    if content.start_with?('---')
      parts = content.split("---", 3)
      frontmatter = parts[1]
      body = parts[2]
      
      # Check if permalink already exists
      return if frontmatter.include?('permalink:')
      
      # Generate permalink from file path
      relative_path = file.sub(@docs_dir + '/', '')
      permalink = generate_permalink(relative_path)
      
      # Add permalink to frontmatter (before closing ---)
      new_frontmatter = frontmatter.rstrip + "\npermalink: #{permalink}\n"
      new_content = "---#{new_frontmatter}---#{body}"
      
      # Write back to file
      File.write(file, new_content)
      puts "Added permalink to #{relative_path}: #{permalink}"
    end
  end

  def generate_permalink(relative_path)
    # Remove .md extension
    path_without_ext = relative_path.sub(/\.md$/, '')
    
    # Ensure it starts and ends with /
    permalink = "/#{path_without_ext}/"
    permalink = permalink.gsub(/\/+/, '/') # Remove double slashes
    permalink
  end
end

# Run the generator
docs_dir = File.expand_path(File.join(__dir__, '..'))
generator = PermalinkGenerator.new(docs_dir)
generator.generate

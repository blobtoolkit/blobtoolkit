#!/usr/bin/env ruby

require 'fileutils'

# Automatically add missing titles to markdown files that already have frontmatter
class FrontmatterGenerator
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
    
    # Only process files that already have frontmatter (start with ---)
    # This ensures we don't add frontmatter to files not meant to be Jekyll pages
    unless content.start_with?('---')
      return
    end
    
    # Parse frontmatter
    parts = content.split("---", 3)
    frontmatter = parts[1]
    body = parts[2]
    
    # Check if title already exists
    return if frontmatter.include?('title:')
    
    # Extract title from first H1 header in the body
    title = extract_title_from_body(body)
    return unless title # Skip if no H1 header found
    
    # Add title to frontmatter (before closing ---)
    new_frontmatter = frontmatter.rstrip + "\ntitle: #{title}\n"
    new_content = "---#{new_frontmatter}---#{body}"
    
    # Write back to file
    File.write(file, new_content)
    relative_path = file.sub(@docs_dir + '/', '')
    puts "Added title to #{relative_path}: #{title}"
  end

  def extract_title_from_body(body)
    # Find the first H1 header
    body.each_line do |line|
      if line.start_with?('# ')
        # Extract text after '# ' and before any trailing whitespace
        title = line.sub(/^#\s+/, '').strip
        return title unless title.empty?
      end
    end
    nil
  end
end

# Run the generator
docs_dir = File.expand_path(File.join(__dir__, '..'))
generator = FrontmatterGenerator.new(docs_dir)
generator.generate

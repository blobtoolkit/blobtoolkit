#!/usr/bin/env ruby

require 'yaml'
require 'fileutils'

# Generate navigation structure from directory layout
class NavGenerator
  def initialize(docs_dir)
    @docs_dir = docs_dir
    @data_dir = File.join(@docs_dir, '_data')
    FileUtils.mkdir_p(@data_dir)
  end

  def generate
    nav_data = {
      'main' => generate_main_nav,
      'sections' => generate_section_navs
    }

    output_file = File.join(@data_dir, 'navigation.yml')
    File.write(output_file, YAML.dump(nav_data))
    puts "✓ Navigation generated: #{output_file}"
    nav_data
  end

  private

  def generate_main_nav
    # Top-level navigation from site config or directory structure
    dirs = Dir.entries(@docs_dir)
      .reject { |d| d.start_with?('.') || d.start_with?('_') || %w[assets scripts].include?(d) }
      .select { |d| File.directory?(File.join(@docs_dir, d)) }
      .sort

    dirs.map do |dir|
      title = dir.split('-').map(&:capitalize).join(' ')
      {
        'title' => title,
        'url' => "/#{dir}/",
        'folder' => dir
      }
    end
  end

  def generate_section_navs
    sections = {}
    main_dirs = Dir.entries(@docs_dir)
      .reject { |d| d.start_with?('.') || d.start_with?('_') || %w[assets scripts].include?(d) }
      .select { |d| File.directory?(File.join(@docs_dir, d)) }

    main_dirs.each do |section_dir|
      section_path = File.join(@docs_dir, section_dir)
      section_nav = []

      # Find all markdown files in section (exclude index.md)
      md_files = Dir.glob(File.join(section_path, '*.md'))
        .reject { |f| File.basename(f) == 'index.md' }
        .sort

      # Add flat markdown files first
      md_files.each do |file|
        filename = File.basename(file, '.md')
        title = extract_title(file) || filename.split('-').map(&:capitalize).join(' ')
        
        section_nav << {
          'title' => title,
          'url' => "/#{section_dir}/#{filename}/",
          'active_pattern' => "/#{section_dir}/#{filename}"
        }
      end

      # Find subdirectories
      subdirs = Dir.entries(section_path)
        .reject { |d| d.start_with?('.') }
        .select { |d| File.directory?(File.join(section_path, d)) }
        .sort

      # Process subdirectories as subsections
      subdirs.each do |subdir|
        subdir_path = File.join(section_path, subdir)
        subdir_index = File.join(subdir_path, 'index.md')
        
        # Only treat as subsection if it has an index.md
        next unless File.exist?(subdir_index)
        
        subdir_title = extract_title(subdir_index) || subdir.split('-').map(&:capitalize).join(' ')
        
        subsection = {
          'title' => subdir_title,
          'url' => "/#{section_dir}/#{subdir}/",
          'active_pattern' => "/#{section_dir}/#{subdir}"
        }

        # Find pages in this subsection (exclude index.md)
        subdir_md_files = Dir.glob(File.join(subdir_path, '*.md'))
          .reject { |f| File.basename(f) == 'index.md' }
          .sort

        # Only add children array if there are actual child pages
        if subdir_md_files.any?
          subsection['children'] = []
          
          subdir_md_files.each do |file|
            filename = File.basename(file, '.md')
            title = extract_title(file) || filename.split('-').map(&:capitalize).join(' ')
            
            subsection['children'] << {
              'title' => title,
              'url' => "/#{section_dir}/#{subdir}/#{filename}/",
              'active_pattern' => "/#{section_dir}/#{subdir}/#{filename}"
            }
          end
        end

        section_nav << subsection
      end

      sections[section_dir] = section_nav
    end

    sections
  end

  def extract_title(file_path)
    return nil unless File.exist?(file_path)
    
    content = File.read(file_path)
    # Look for front matter title
    if content.match?(/^---/)
      match = content.match(/^---\n(.*?)\n---/m)
      return nil unless match
      
      front_matter = match[1]
      title_match = front_matter.match(/^title:\s*(.+)$/i)
      return title_match[1].strip if title_match
    end
    
    nil
  end
end

# Run the generator
if __FILE__ == $PROGRAM_NAME
  docs_dir = File.expand_path('../', __dir__)
  generator = NavGenerator.new(docs_dir)
  generator.generate
end

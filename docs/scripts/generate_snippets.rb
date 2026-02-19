#!/usr/bin/env ruby

require 'json'
require 'nokogiri'
require 'fileutils'

# Check if _site exists
unless Dir.exist?('_site')
  puts "⚠️  _site directory not found. Run Jekyll build first."
  exit 0
end

# Get all HTML files from _site
snippets = []

Dir.glob('_site/**/*.html').each do |file|
  # Skip 404 and feed files
  next if file.include?('404.html') || file.include?('feed.xml')
  
  # Read and parse the HTML
  html_content = File.read(file)
  doc = Nokogiri::HTML(html_content)
  
  # Extract text content from the main site-content div
  content_div = doc.css('.site-content').first
  next unless content_div
  
  # Remove all heading elements
  content_div.css('h1, h2, h3, h4, h5, h6').each(&:remove)
  
  # Get all paragraph elements
  paragraphs = content_div.css('p')
  
  snippet = ''
  if paragraphs.any?
    # Get first paragraph text
    first_para = paragraphs.first.text.strip
    
    # Extract first sentence (ending with . ! or ?)
    if first_para =~ /^([^.!?]+[.!?])/
      snippet = $1.strip
    else
      # If no sentence ending, take first 150 chars
      snippet = first_para.length > 150 ? first_para[0..149] + '...' : first_para
    end
  else
    # Fallback: get cleaned text
    text = content_div.text.gsub(/\s+/, ' ').strip
    snippet = text.length > 150 ? text[0..149] + '...' : text
  end
  
  # Convert file path to URL
  url = '/' + file.sub('_site/', '').sub('index.html', '').sub('.html', '/')
  url = '/' if url == '/'
  
  snippets << {
    'id' => url,
    'snippet' => snippet
  }
end

# Write JSON file
FileUtils.mkdir_p('_site')
output = { 'snippets' => snippets }
File.write('_site/snippets.json', JSON.pretty_generate(output))

puts "✅ Generated #{snippets.length} snippets"

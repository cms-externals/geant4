#!/usr/bin/env python3

import re
import sys
import os

def convert_section_line(text, level):
    """
    Helper function to convert the section line
    """
    words = re.split(r'[\s]+', text.strip())    
    # print(words)    
    new_text = ''
    for i in range(0, len(words)):
        if i == 0:
            # replace the section/subsection keyword with ##...
            for j in range(0,level):
                new_text += '#'   
        if i > 1:
            new_text += ' '
        if i > 1:
            new_text += words[i]
    return new_text
    
def convert_doxygen_to_markdown(doxygen_text):
    """
    Converts a string with Doxygen-like formatting to Markdown.

    Supports:
    - \\section {id} {title} -> ## title (id skipped, handles multi-word IDs heuristically)
    - \\subsection {id} {title} -> ### title (id skipped, handles multi-word IDs heuristically)
    - \\subsubsection {id} {title} -> #### title (id skipped, handles multi-word IDs heuristically)
    - \\verbatim ... \\endverbatim -> ```\n...\n```
    - ///\\file -> (line is skipped entirely)
    - ///\\brief -> (line is skipped entirely if it follows ///\\file)
    - /*! \\page -> \\page (line is kept, only prefix removed)
    - Removes lines containing '*/' (possibly with leading/trailing whitespace)
    """
    markdown_text = []
    lines = doxygen_text.splitlines()
    in_verbatim_block = False
    prev_line_was_file_directive = False

    # --- New Logic: Find the starting point (line with \page) ---
    start_processing_index = 0
    found_page_directive = False
    for i, line in enumerate(lines):
        # Check for both '/*! \page' and '\page' as the target for the starting point
        if r'\page' in line or r'/*! \page' in line:
            start_processing_index = i
            found_page_directive = True
            break

    # If '\page' was not found, process all lines from the beginning.
    # Otherwise, start from the line where '\page' was found.
    lines_to_process = lines[start_processing_index:]

    for line in lines_to_process:
        current_line_stripped = line.strip()

        # Apply the /*! \page replacement first, as this line might be the starting point
        line = line.replace(r'/*! \page', r'\page')
        current_line_stripped = line.strip() # Re-strip after potential replacement

        # Priority 1: Handle lines that should always be skipped or affect state for the *next* line
        # Check if the line (after stripping whitespace) contains '*/'
        if '*/' in current_line_stripped:
            prev_line_was_file_directive = False # Reset flag if '*/' interrupts sequence
            continue # Skip this line entirely

        # Check if the line starts with '///\file'
        if current_line_stripped.startswith(r'///\file'):
            prev_line_was_file_directive = True # Set flag for the *next* line
            continue # Skip this line entirely

        # Priority 2: Handle ///\brief conditionally
        # This check must come AFTER the ///\file check, as it depends on prev_line_was_file_directive
        if prev_line_was_file_directive and current_line_stripped.startswith(r'///\brief'):
            # This line is skipped because the condition is met.
            # Reset the flag immediately as this brief has been handled.
            prev_line_was_file_directive = False
            continue # Skip this line entirely

        # Priority 3: For any other line that was not skipped by the above conditions,
        # ensure the 'prev_line_was_file_directive' flag is reset.
        # This ensures that a ///\brief line is only skipped if it's *immediately* after a ///\file.
        # If any other line type comes between ///\file and ///\brief, the brief should NOT be skipped.
        prev_line_was_file_directive = False # Reset for generic lines

        # Remaining processing logic
        if in_verbatim_block:
            if current_line_stripped == r'\endverbatim':
                markdown_text.append('```')
                in_verbatim_block = False
            else:
                markdown_text.append(line)
        else:
            if current_line_stripped.startswith(r'\section'):
                markdown_text.append(convert_section_line(current_line_stripped, 2))
            elif current_line_stripped.startswith(r'\subsection'):
                # Similar logic for subsection
                markdown_text.append(convert_section_line(current_line_stripped, 3))
            elif current_line_stripped.startswith(r'\subsubsection'):
                # Similar logic for subsection
                markdown_text.append(convert_section_line(current_line_stripped, 4))
            elif current_line_stripped == r'\verbatim':
                markdown_text.append('```')
                in_verbatim_block = True
            else:
                markdown_text.append(line)

    return "\n".join(markdown_text)

def batch_convert_doxygen_readmes(start_directory):
    """
    Searches for all files named *.README.txt in the given directory and its subdirectories,
    converts them to Markdown, and saves them as *.README.md in the same directory.
    """
    print(f"Starting batch conversion in directory: {start_directory}")
    converted_count = 0
    skipped_count = 0

    for root, _, files in os.walk(start_directory):
        for filename in files:
            if filename.endswith(".README.txt"):
                input_filepath = os.path.join(root, filename)
                output_filename = filename.replace(".README.txt", "README.md")
                output_filepath = os.path.join(root, output_filename)

                print(f"Processing: {input_filepath}")
                try:
                    with open(input_filepath, 'r', encoding='utf-8') as infile:
                        doxygen_content = infile.read()

                    markdown_output = convert_doxygen_to_markdown(doxygen_content)

                    with open(output_filepath, 'w', encoding='utf-8') as outfile:
                        outfile.write(markdown_output)
                    print(f"  Converted to: {output_filepath}")
                    converted_count += 1
                except Exception as e:
                    print(f"  Error converting '{input_filepath}': {e}")
                    skipped_count += 1
    print(f"\nBatch conversion complete. Converted {converted_count} files, skipped {skipped_count} files due to errors.")


# --- Example Usage ---
if __name__ == "__main__":
    # Check for batch conversion mode first, as it also has 3 arguments
    if len(sys.argv) == 3 and sys.argv[1] == "--batch":
        start_directory = sys.argv[2]
        if not os.path.isdir(start_directory):
            print(f"Error: Directory '{start_directory}' not found.")
            sys.exit(1)
        batch_convert_doxygen_readmes(start_directory)
    # Then check for single file conversion mode
    elif len(sys.argv) == 3:
        input_filename = sys.argv[1]
        output_filename = sys.argv[2]

        try:
            with open(input_filename, 'r', encoding='utf-8') as infile:
                doxygen_content = infile.read()
        except FileNotFoundError:
            print(f"Error: Input file '{input_filename}' not found.")
            sys.exit(1)
        except Exception as e:
            print(f"Error reading input file: {e}")
            sys.exit(1)

        markdown_output = convert_doxygen_to_markdown(doxygen_content)

        try:
            with open(output_filename, 'w', encoding='utf-8') as outfile:
                outfile.write(markdown_output)
            print(f"Successfully converted '{input_filename}' to '{output_filename}'.")
        except Exception as e:
            print(f"Error writing to output file: {e}")
            sys.exit(1)
    # If neither of the above, print usage
    else:
        print("Usage for single file: python doxygen_converter.py <input_file.dox> <output_file.md>")
        print("Usage for batch conversion: python doxygen_converter.py --batch <start_directory>")
        print("Example (single): python doxygen_converter.py input.dox output.md")
        print("Example (batch): python doxygen_converter.py --batch ./my_doxygen_docs")
        sys.exit(1)

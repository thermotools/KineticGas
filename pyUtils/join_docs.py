"""
Scripts to generate readme's and wiki-documentation files by joining various markdown files in KineticGas/docs/markdown

Intent: Several "submodules" of documentation are used both in the wiki, the main github README, the pypi README and
        possibly other places. Each of these "submodules" should be contained in a single markdown file in the
        thermopack/doc/markdown directory. Each function in this file generates a single markdown file by joining
        the appropriate submodules, and prepending the header generated by the get_header(files) function.

Usage: To add new documentation, create a new markdown file in thermopoack/doc/markdown, and add the filename
        (sans the file ending) to the appropriate `files` lists in the functions in this file.
"""
import os
from datetime import datetime
from tools import write_file

MARKDOWN_DIR = os.path.dirname(__file__) + '/../docs/'
KINETICGAS_ROOT = os.path.dirname(__file__) + '/..'

def print_finished_report(header, out_file_path):
    printcolwidth = 100
    print('-' * printcolwidth)
    print('Wrote documentation file with the header :')
    print('#' * printcolwidth)
    print(header)
    print('#' * printcolwidth)
    print('-' * printcolwidth)
    print('To:', out_file_path)
    print('-' * printcolwidth)
    print()

def gen_file_str(files):
    out_file_str = ''
    for file in files:
        file_path = MARKDOWN_DIR + file + '.md'

        with open(file_path, 'r') as in_file:
            out_file_str += in_file.read() + '\n\n'

    return out_file_str

def get_header(files):
    header = '<!--- \n'
    header += 'Generated at: ' + datetime.today().isoformat() + '\n'
    header += 'This is an auto-generated file, generated using the script at KineticGas/docs/join_docs.py\n'
    header += 'The file is created by joining the contents of the files\n'
    header += '    KineticGas/docs/markdown/\n'
    for fname in files:
        header += '    ' * 2 + fname + '.md\n'
    header += '--->\n'
    return header

def write_pypi_readme():
    files = ['readme_parts/header', 'readme_parts/toc_pypi', 'metapages/cite_acknowl_licence', 'readme_parts/pypi_structure',
             'vCurrent/getting_started_py', 'vCurrent/fluid_identifiers']
    header = get_header(files)

    out_file_str = header + gen_file_str(files)
    out_file_path = KINETICGAS_ROOT + '/pykingas/README.md'
    write_file(out_file_path, out_file_str)

def write_github_readme():
    files = ['readme_parts/header', 'readme_parts/toc_github', 'metapages/cite_acknowl_licence', 'vCurrent/source_build',
            'vCurrent/getting_started_py', 'vCurrent/getting_started_cpp', 'vCurrent/advanced', 'vCurrent/structure',
             'vCurrent/fluid_identifiers']
    header = get_header(files)

    out_file_str = header + gen_file_str(files)
    out_file_path = KINETICGAS_ROOT + '/README.md'
    write_file(out_file_path, out_file_str)

if __name__ == '__main__':
    write_pypi_readme()
    write_github_readme()


'''
Functions to auto-generate markdown documentation by parsing docstrings

Intent : As long as we keep docstrings consistently formatted, this ensures that we can quickly generate up-to-date
        documentation that is easily accessible. It also ensures that documentation only needs to be updated one place,
        and that place is in the codebase itself.

Usage : Current functionality is designed to parse the docstrings of a given class. It assumes that the docstrings of
        the methods in the class are formatted as

        def myfunc(self, p1, p2, p3, p4=None, p5=<something>, ...):
            """Section Name
            Description of what this function does (a kind of header). We can write lots of stuff here
            NOTE the double lineshift here

            Args:
                 p1 (int) : The lineshift before 'Args' is necessary.
                 p2 (float) : The colons here are also necessary.
                 p3 (bool) : SomethingSomething
                 p4 (Optional, list) : etc.
            Returns:
                (float) : The colon here is also necessary.
            """
-----------------------------------------------------------------------------------------------------------------------

NOTE: The procedure described below is the general procedure for generating documentation for large classes
        with many methods doing different things. For simple classes, consider using the function
        `basic_class_to_markdown` to do most of the job.

-----------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------
Now for the generic procedure:

    1 :  Use the function get_autogen_header(<classname>) to generate a comment with a timestamp and a description of
         how the markdown file was generated (using this module)

    2 : Define a list of sections that the methods in the class should be divided into, see : thermo_to_markdown() for
        an example. The sections should equate to the 'Section Name' in the above formatting example. The sections in
        the markdown file will be ordered according to the order of this list.

        NOTE: The sections "Other" and "Deprecated" are automatically generated.
            If the section name in a methods docstring does not match any of the supplied sections, it is added to the
            "Other" section, and a warning is issued
            If a method has the Section Name "Deprecated" it is added to the automatically generated "Deprecated" section.
            To add other automatically generated sections, add them in the function `get_automatic_sections`

        NOTE: A method can be added to several sections, by separating the "Section Name"'s with an ampersand as
            def my_func(self):
                """Section 1 & Section 5 & Section 6 & ...
                This is documentation for a method that belongs to many sections
                """

    3 : Define a dict of section name => section header mappings. This maps the section names in the docstrings to
        the section headers in the markdown file.

    4 : Define a dict of section name => section intro mappings. This maps the section names to an intro text that is
        displayed at the start of the section.

    5 : Use
            my_methods = inspect.getmembers(<MyClass>, predicate=inspect.isfunction)
        to extract the list of class methods.

    6 : Use
            method_dict = split_methods_by_section(sections, my_methods)
        to split the methods into a single list for each method, placed in a dict

    7 : Use
            toc_markdown_text = get_toc(sections, section_headers, method_dict)
        to generate a string, consisting of the markdown table of contents (with links)

    8 : Use
            contents = get_markdown_contents(sections, section_headers, section_intro, method_dict)
        to generate the contents of the file

    9 : Concatenate the strings you have generated and use `write_file` to write to file.
    Note: USE The method `write_file`: It checks for content changes and does not update the file if no
        other content than the timestamp has changed. This prevents you from needing to push a bunch of changes that
        are only timestamp changes, and makes our history much more clean.
'''
import copy
import inspect
from warnings import warn
from datetime import datetime
from tools import remove_illegal_link_chars, check_is_changed, write_file, KINETICGAS_ROOT, MARKDOWN_DIR
from pykingas.py_KineticGas import py_KineticGas
from pykingas.MieType import MieType
from pykingas.MieKinGas import MieKinGas
from pykingas.HardSphere import HardSphere


def get_autogen_header(classname):
    header = f'''---
layout: default
version: 
title: Methods in the {classname} class
permalink: /vcurrent/{classname}_methods.html
---\n\n'''
    header += '<!--- \n'
    header += 'Generated at: ' + datetime.today().isoformat() + '\n'
    header += 'This is an auto-generated file, generated using the script at ' \
              'KineticGas/pyUtils/markdown_from_docstrings.py\n'
    header += 'The file is created by parsing the docstrings of the methods in the \n'
    header += classname + ' class. For instructions on how to use the parser routines, see the\n' \
                          'file KineticGas/pyUtils/markdown_from_docstrings.py'
    header += '--->\n\n'
    return header


def to_markdown(methods):
    """
    Generate markdown text from the docstring of the methods in the list methods

    Args:
        methods list<tuple<str, function>> : A list of tuples, with the first element being the name of the method,
                                            and the second element being the method itself.
    Returns:
        (str) : A markdown - formatted string with the contents of the docstring.
    """
    md_text = ''
    for name, meth in methods:
        docparts = meth.__doc__.split('\n\n')
        header_lines = [line.strip() for line in docparts[0].split('\n')]
        header = ' '.join(header_lines[1:]) # Cutting out the section identifier
        header.replace('\n', ' ')

        if len(docparts) > 1:
            content = '\n\n'.join(docparts[1:])
        else:
            content = ''

        content_lines = [line.strip() for line in content.split('\n')]

        md_text += '### `' + name + str(inspect.signature(meth)) + '`\n'
        md_text += header + '\n\n'

        pad = '&nbsp;' * 4 + ' '
        endl = '\n\n'
        for line in content_lines:
            if ('args:' in line.lower()) or ('returns:' in line.lower()) or ('raises:' in line.lower()):
                md_text += '#### ' + line + endl

            elif ':' in line:
                line = line.split(':')
                md_text += pad + '**' + line[0] + ':** ' + endl
                md_text += 2 * pad + ':'.join(line[1:]) + endl
            else:
                md_text += 2 * pad + line + endl

    return md_text

def split_methods_by_section(sections, methods):
    """
    Organise the various methods of a class into sections, determined by the id 'Section Name' as in the example at
    the top of this file. Warns if there are methods for which no matching section was found.
    The section 'Other' is automatically generated if there are methods with an identifier not matching any of the
    identifiers in 'sections'. The section 'Deprecated' is automatically included in the 'sections' list.

    Args:
        sections (list<str>) : The name of each section, corresponding to the first line in each docstring.
        methods (list<tuple<str, function>>) : List of tuples, with the first element corresponding to the method
                                                names, and the second to the methods themself.
    Returns:
        dict : With keys corresponding to the sections, and value being a list of the (name, method) tuples with
                the first line of the docstring matching the section (key).
    """

    if 'deprecated' not in [s.lower() for s in sections]:
        sections.append('Deprecated')

    method_dict = {}
    for name, meth in methods:
        # A method can be added to several sections, by giving it the header
        # def myfunc():
        #   """Section1 & Section2 & Section5 ... """
        method_sections = [m.strip() for m in meth.__doc__.split('\n')[0].lower().split('&')] # extracting the section names as described above
        method_has_section = False
        for sec in sections:
            if sec.lower() in method_sections:
                method_has_section = True
                if sec in method_dict.keys():
                    method_dict[sec].append((name, meth))
                else:
                    method_dict[sec] = [(name, meth)]
        if method_has_section is False: # Generate a section called "Other" and add the method to that section. Warn that this is done.
            warn('Method : ' + name + " did not contain a section on the first line of its docstring! Adding to 'Other'",
                 SyntaxWarning, stacklevel=3)
            warn('Section name in docstring is : '+ meth.__doc__.split('\n')[0].lower(), SyntaxWarning, stacklevel=3)

            if 'Other' in method_dict.keys():
                method_dict['Other'].append((name, meth))
            else:
                method_dict['Other'] = [(name, meth)]

        # print('other is : ', method_dict['Other'] if 'Other' in method_dict.keys() else None)
    return method_dict

def get_automatic_sections(sections, section_headers, section_intro, method_dict):
    sections = copy.deepcopy(sections)
    section_headers = copy.deepcopy(section_headers)
    section_intro = copy.deepcopy(section_intro)

    if 'Other' in method_dict.keys():
        if 'Other' not in sections:
            sections.append('Other')
        if 'Other' not in section_headers.keys():
            section_headers['Other'] = 'Other'
        if 'Other' not in section_intro.keys():
            section_intro['Other'] = 'Methods that do not have a section identifier in their docstring.'

    if 'Deprecated' in method_dict.keys():
        if 'Deprecated' not in sections:
            sections.append('Deprecated')
        if 'Deprecated' not in section_headers.keys():
            section_headers['Deprecated'] = 'Deprecated methods'
        if 'Deprecated' not in section_intro.keys():
            section_intro['Deprecated'] = 'Deprecated methods are not maintained, and may be removed in the future.'

    return sections, section_headers, section_intro

def get_toc(sections, section_headers, method_dict, is_subsection=False):
    """
    Generate a table of contents with links to the sections, using the names in section_headers.
    Note: Uses the function tools.remove_illegal_link_chars() to sanitize raw strings, such that they can be used
    in markdown links.

    Args:
        sections (list<str>) : List of the section names, corresponding to the first line in the docstrings
        section_headers (dict<str, str>) : Dict mapping the section names to the section headers in the markdown file
        method_dict (dict<str, list<tuple<str, function>>>) : Mapping the section names to the list of (method name,
                                                                function) tuples that are in the section.
        is_subsection (bool, optional) : Whether to include the automatically generated sections in the toc.
                                        Should be set to False when generating the toc for the whole file, and set to
                                        True when generating toc for subsections.
    Returns:
        str : The string representation of the table of contents to be written to the markdown file.
    """

    if is_subsection is False:
        sections, section_headers, _ = get_automatic_sections(sections, section_headers, {}, method_dict)

    toc_text = '## Table of contents\n'
    for sec in sections:
        if sec not in method_dict.keys():
            continue

        sec_name = section_headers[sec]
        sec_id = remove_illegal_link_chars(sec_name)
        toc_text += '  * [' + sec_name + '](#' + sec_id + ')\n'

        for meth in method_dict[sec]:
            method_name = meth[0].replace('__', '\_\_')
            method_id = meth[0] + str(inspect.signature(meth[1]))
            method_id = remove_illegal_link_chars(method_id)
            toc_text += '    * [' + method_name + '](#' + method_id + ')\n'

    return toc_text + '\n'

def get_markdown_contents(sections, section_headers, section_intro, method_dict):
    """
    Iterate through the sections, generate the markdown documentation for all methods in method_dict, and join them
    while adding section headers and section introductions

    Args:
        sections (list<str>) : List of the section names, corresponding to the first line in the docstrings
        section_headers (dict<str, str>) : Dict mapping the section names to the section headers in the markdown file
        section_intro (dict<str, str>) : Dict mapping the section names to the section introduction text.
        method_dict (dict<str, list<tuple<str, function>>>) : Mapping the section names to the list of (method name,
                                                                function) tuples that are in the section.
    Returns:
        str : The markdown text corresponding to the main contents of the file.
    """

    sections, section_headers, section_intro = get_automatic_sections(sections, section_headers, section_intro, method_dict)

    md_text = ''
    for sec in sections:
        if sec not in method_dict.keys():
            continue
        md_text += '## ' + section_headers[sec] + '\n\n'
        md_text += section_intro[sec] + '\n\n'
        md_text +=  '#' + get_toc([sec], section_headers, method_dict, is_subsection=True) + '\n'
        md_text += to_markdown(method_dict[sec])

    return md_text

def basic_class_to_markdown(classname, eosname, methods, intro_text=None, inherits=None):
    """
    Generate markdown documentation file for a class that implements only Constructor and unility methods.
    """

    sections = ['Constructor',
                'Utility']

    section_headers = {'Constructor': 'Constructor',
                       'Utility': 'Utility methods'}

    section_intro = {'Constructor': f'Methods to initialise {eosname} model.',
                     'Utility': 'Set- and get methods for interaction parameters, mixing parameters ...'}

    method_dict = split_methods_by_section(sections, methods)

    ofile_text = get_autogen_header(classname)
    if intro_text is None:
        if inherits is None:
            ofile_text += f'The `{classname}` class, found in `pykingas/{classname}.py`, is the interface to the \n' \
                          f'{eosname} Model. This class implements utility methods to access parameters etc.\n\n'
        else:
            ofile_text += f'The `{classname}` class, found in `pykingas/{classname}.py`, inherrits ' \
                          f'from the {inherits} class, and  is the interface to the \n' \
                          f'{eosname} Model. This class implements utility methods to access mixing parameters etc.\n\n'
    else:
        ofile_text += intro_text
    ofile_text += get_toc(sections, section_headers, method_dict)
    ofile_text += get_markdown_contents(sections, section_headers, section_intro, method_dict)

    filename = f'{classname}_methods.md'
    write_file(MARKDOWN_DIR + filename, ofile_text)

def pyKineticGas_to_markdown():
    """
    Generate markdown documentation file for the thermo class.
    """

    sections = ['Constructor',
                'TV-property',
                'Tp-property',
                'FOR-Transform',
                'cpp-interface',
                'Utility']

    section_headers = {'Constructor' : 'The constructor',
                        'Utility' : 'Utility methods',
                       'TV-property' : 'TV-property interfaces',
                       'Tp-property' : 'Tp-property interfaces',
                       'FOR-Transform' : 'Frame of Reference transformations',
                       'cpp-interface' : 'Interfaces to C++ methods'}

    section_intro = {'Constructor' : 'The constructor',
                    'Utility': 'Methods for various helper computations',
                       'TV-property': 'Computing properties as a function of temperature and volume.',
                       'Tp-property': 'Computing properties as a function of temperature and pressure. Simply forwards '
                                      'calls to the appropriate TV-property methods after computing molar volume using '
                                      'the internal equation of state object.',
                       'FOR-Transform': 'Generate matrices for Frame of Reference transformations. See the supporting'
                                        'material of [Revised Enskog Theory for Mie fluids](https://doi.org/10.1063/5.0149865) '
                                        'for details.',
                       'cpp-interface': 'Lightweight wrappers for the most commonly used C++ side methods.'}

    kingas_methods = inspect.getmembers(py_KineticGas, predicate=inspect.isfunction)
    method_dict = split_methods_by_section(sections, kingas_methods)

    ofile_text = get_autogen_header('py_KineticGas')
    ofile_text += 'The `py_KineticGas` class, found in `pykingas/py_KineticGas.py`, is the core of the KineticGas Python interface. ' \
                  'All Revised Enskog Theory models inherit from `py_KineticGas`. This is the class that contains the interface to all ' \
                  'practical calculations that can be done from the python-side of KineticGas. Derived classes only implement ' \
                  'specific functions for parameter handling etc.\n\n'
    ofile_text += get_toc(sections, section_headers, method_dict)
    ofile_text += get_markdown_contents(sections, section_headers, section_intro, method_dict)

    filename = 'py_KineticGas_methods.md'
    write_file(MARKDOWN_DIR + filename, ofile_text)


def mietype_to_markdown():
    classname = 'MieType'
    eosname = 'Mie-Type'

    class_methods = inspect.getmembers(MieType, predicate=inspect.isfunction)
    parent_methods = inspect.getmembers(py_KineticGas, predicate=inspect.isfunction)
    specific_methods = sorted(list(set(class_methods) - set(parent_methods)))

    basic_class_to_markdown(classname, eosname, specific_methods, inherits='py_KineticGas')

def miekingas_to_markdown():
    classname = 'MieKinGas'
    eosname = 'RET-Mie'

    class_methods = inspect.getmembers(MieKinGas, predicate=inspect.isfunction)
    parent_methods = inspect.getmembers(MieType, predicate=inspect.isfunction)
    specific_methods = sorted(list(set(class_methods) - set(parent_methods)))

    basic_class_to_markdown(classname, eosname, specific_methods, inherits='MieType')

def hardsphere_to_markdown():
    classname = 'HardSphere'
    eosname = 'RET-HS'

    class_methods = inspect.getmembers(HardSphere, predicate=inspect.isfunction)
    parent_methods = inspect.getmembers(py_KineticGas, predicate=inspect.isfunction)
    specific_methods = sorted(list(set(class_methods) - set(parent_methods)))

    basic_class_to_markdown(classname, eosname, specific_methods, inherits='py_KineticGas')


if __name__ == '__main__':
    pyKineticGas_to_markdown()
    mietype_to_markdown()
    miekingas_to_markdown()
    hardsphere_to_markdown()
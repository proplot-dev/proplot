# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This directive will produce an "autosummary"-style table for public
attributes of a specified module. See the `sphinx.ext.autosummary`_ extension
for details on this process. The main difference from the `autosummary`_
directive is that `autosummary`_ requires manually inputting all attributes
that appear in the table, while this captures the entries automatically.

This directive requires a single argument that must be a module or
package.

It also accepts any options supported by the `autosummary`_ directive-
see `sphinx.ext.autosummary`_ for details. It also accepts some additional
options:

    * ``:classes-only:``
        If present, the autosummary table will only contain entries for
        classes. This cannot be used at the same time with
        ``:functions-only:`` or ``:variables-only:``.

    * ``:functions-only:``
        If present, the autosummary table will only contain entries for
        functions. This cannot be used at the same time with
        ``:classes-only:`` or ``:variables-only:``.

    * ``:variables-only:``
        If present, the autosummary table will only contain entries for
        variables (everything except functions and classes). This cannot
        be used at the same time with ``:classes-only:`` or
        ``:functions-only:``.

    * ``:skip: obj1, [obj2, obj3, ...]``
        If present, specifies that the listed objects should be skipped
        and not have their documentation generated, nor be included in
        the summary table.

    * ``:allowed-package-names: pkgormod1, [pkgormod2, pkgormod3, ...]``
        Specifies the packages that functions/classes documented here are
        allowed to be from, as comma-separated list of package names. If not
        given, only objects that are actually in a subpackage of the package
        currently being documented are included.

    * ``:inherited-members:`` or ``:no-inherited-members:``
        The global sphinx configuration option ``automodsumm_inherited_members``
        decides if members that a class inherits from a base class are included
        in the generated documentation. The flags ``:inherited-members:`` or
        ``:no-inherited-members:`` allows overrriding this global setting.

This extension also adds two sphinx configuration options:

* ``automodsumm_writereprocessed``
    Should be a bool, and if ``True``, will cause `automodsumm`_ to write files
    with any ``automodsumm`` sections replaced with the content Sphinx
    processes after ``automodsumm`` has run.  The output files are not
    actually used by sphinx, so this option is only for figuring out the
    cause of sphinx warnings or other debugging.  Defaults to ``False``.

* ``automodsumm_inherited_members``
    Should be a bool and if ``True``, will cause `automodsumm`_ to document
    class members that are inherited from a base class. This value can be
    overriden for any particular automodsumm directive by including the
    ``:inherited-members:`` or ``:no-inherited-members:`` options.  Defaults to
    ``False``.

.. _sphinx.ext.autosummary: http://sphinx-doc.org/latest/ext/autosummary.html
.. _autosummary: http://sphinx-doc.org/latest/ext/autosummary.html#directive-autosummary

.. _automod-diagram:

automod-diagram directive
=========================

This directive will produce an inheritance diagram like that of the
`sphinx.ext.inheritance_diagram`_ extension.

This directive requires a single argument that must be a module or
package. It accepts no options.

.. note::
    Like 'inheritance-diagram', 'automod-diagram' requires
    `graphviz <http://www.graphviz.org/>`_ to generate the inheritance diagram.

.. _sphinx.ext.inheritance_diagram: http://sphinx-doc.org/latest/ext/inheritance.html
"""

import abc
import inspect
import os
import re
import io
from distutils.version import LooseVersion

from sphinx import __version__
from sphinx.ext.autosummary import Autosummary
from sphinx.ext.inheritance_diagram import InheritanceDiagram
from docutils.parsers.rst.directives import flag

from .utils import find_mod_objs, cleanup_whitespace

__all__ = ['Automoddiagram', 'Automodsumm', 'automodsumm_to_autosummary_lines',
           'generate_automodsumm_docs', 'process_automodsumm_generation']

SPHINX_LT_16 = LooseVersion(__version__) < LooseVersion('1.6')
SPHINX_LT_17 = LooseVersion(__version__) < LooseVersion('1.7')

api_class_methods = ['__init__', '__call__',
        '__getitem__', '__setitem__', '__setattr__',
        '__getattribute__']
api_ignore_methods = []

def _str_list_converter(argument):
    """
    A directive option conversion function that converts the option into a list
    of strings. Used for 'skip' option.
    """
    if argument is None:
        return []
    else:
        return [s.strip() for s in argument.split(',')]


class Automodsumm(Autosummary):
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = False
    has_content = False
    option_spec = dict(Autosummary.option_spec)
    option_spec['functions-only'] = flag
    option_spec['classes-only'] = flag
    option_spec['variables-only'] = flag
    option_spec['skip'] = _str_list_converter
    option_spec['allowed-package-names'] = _str_list_converter
    option_spec['inherited-members'] = flag
    option_spec['no-inherited-members'] = flag

    def run(self):
        env = self.state.document.settings.env
        modname = self.arguments[0]

        try:
            self.warnings[:] = []
        except AttributeError:
            self.warnings = []

        nodelist = []

        try:
            localnames, fqns, objs = find_mod_objs(modname)
        except ImportError:
            self.warn("Couldn't import module " + modname)
            return self.warnings

        try:
            # set self.content to trick the autosummary internals.
            # Be sure to respect functions-only and classes-only.
            funconly = 'functions-only' in self.options
            clsonly = 'classes-only' in self.options
            varonly = 'variables-only' in self.options
            if [clsonly, funconly, varonly].count(True) > 1:
                self.warning('more than one of functions-only, classes-only, '
                             'or variables-only defined. Ignoring.')
                clsonly = funconly = varonly = False

            skipnames = []
            if 'skip' in self.options:
                option_skipnames = set(self.options['skip'])
                for lnm in localnames:
                    if lnm in option_skipnames:
                        option_skipnames.remove(lnm)
                        skipnames.append(lnm)
                if len(option_skipnames) > 0:
                    self.warn('Tried to skip objects {objs} in module {mod}, '
                              'but they were not present.  Ignoring.'
                              .format(objs=option_skipnames, mod=modname))

            if funconly:
                cont = []
                for nm, obj in zip(localnames, objs):
                    if nm not in skipnames and inspect.isroutine(obj):
                        cont.append(nm)
            elif clsonly:
                cont = []
                for nm, obj in zip(localnames, objs):
                    if nm not in skipnames and inspect.isclass(obj):
                        cont.append(nm)
            elif varonly:
                cont = []
                for nm, obj in zip(localnames, objs):
                    if nm not in skipnames and not (inspect.isclass(obj) or
                                                    inspect.isroutine(obj)):
                        cont.append(nm)
            else:
                cont = [nm for nm in localnames if nm not in skipnames]

            self.content = cont

            # for some reason, even though ``currentmodule`` is substituted in,
            # sphinx doesn't necessarily recognize this fact.  So we just force
            # it internally, and that seems to fix things
            env.temp_data['py:module'] = modname
            env.ref_context['py:module'] = modname

            # can't use super because Sphinx/docutils has trouble return
            # super(Autosummary,self).run()
            nodelist.extend(Autosummary.run(self))

            return self.warnings + nodelist
        finally:  # has_content = False for the Automodsumm
            self.content = []

    def get_items(self, names):
        self.genopt['imported-members'] = True
        return Autosummary.get_items(self, names)


# <-------------------automod-diagram stuff----------------------------------->
class Automoddiagram(InheritanceDiagram):

    option_spec = dict(InheritanceDiagram.option_spec)
    option_spec['allowed-package-names'] = _str_list_converter
    option_spec['skip'] = _str_list_converter

    def run(self):
        try:
            ols = self.options.get('allowed-package-names', [])
            ols = True if len(ols) == 0 else ols  # if none are given, assume only local

            nms, objs = find_mod_objs(self.arguments[0], onlylocals=ols)[1:]
        except ImportError:
            self.warnings = []
            self.warn("Couldn't import module " + self.arguments[0])
            return self.warnings

        # Check if some classes should be skipped
        skip = self.options.get('skip', [])

        clsnms = []
        for n, o in zip(nms, objs):

            if n.split('.')[-1] in skip:
                continue

            if inspect.isclass(o):
                clsnms.append(n)

        oldargs = self.arguments
        try:
            if len(clsnms) > 0:
                self.arguments = [' '.join(clsnms)]
            return InheritanceDiagram.run(self)
        finally:
            self.arguments = oldargs


# <---------------------automodsumm generation stuff-------------------------->
def process_automodsumm_generation(app):
    env = app.builder.env

    filestosearch = []
    for docname in env.found_docs:
        filename = env.doc2path(docname)
        if os.path.isfile(filename):
            filestosearch.append(docname + os.path.splitext(filename)[1])

    liness = []
    for sfn in filestosearch:
        lines = automodsumm_to_autosummary_lines(sfn, app)
        liness.append(lines)
        if app.config.automodsumm_writereprocessed:
            if lines:  # empty list means no automodsumm entry is in the file
                outfn = os.path.join(app.srcdir, sfn) + '.automodsumm'
                with open(outfn, 'w') as f:
                    for l in lines:
                        f.write(l)
                        f.write('\n')

    for sfn, lines in zip(filestosearch, liness):
        suffix = os.path.splitext(sfn)[1]
        if len(lines) > 0:
            generate_automodsumm_docs(
                lines, sfn, app=app, builder=app.builder,
                suffix=suffix, base_path=app.srcdir,
                inherited_members=app.config.automodsumm_inherited_members)

# _automodsummrex = re.compile(r'^(\s*)\.\. automodsumm::\s*([A-Za-z0-9_.]+)\s*'
#                              r'\n\1(\s*)(\S|$)', re.MULTILINE)
_lineendrex = r'(?:\n|$)'
_hdrex = r'^\n?(\s*)\.\. automodsumm::\s*(\S+)\s*' + _lineendrex
_oprex1 = r'(?:\1(\s+)\S.*' + _lineendrex + ')'
_oprex2 = r'(?:\1\4\S.*' + _lineendrex + ')'
_automodsummrex = re.compile(_hdrex + '(' + _oprex1 + '?' + _oprex2 + '*)',
                             re.MULTILINE)


def automodsumm_to_autosummary_lines(fn, app):
    """
    Generates lines from a file with an "automodsumm" entry suitable for
    feeding into "autosummary".

    Searches the provided file for `automodsumm` directives and returns
    a list of lines specifying the `autosummary` commands for the modules
    requested. This does *not* return the whole file contents - just an
    autosummary section in place of any :automodsumm: entries. Note that
    any options given for `automodsumm` are also included in the
    generated `autosummary` section.

    Parameters
    ----------
    fn : str
        The name of the file to search for `automodsumm` entries.
    app : sphinx.application.Application
        The sphinx Application object

    Returns
    -------
    lines : list of str
        Lines for all `automodsumm` entries with the entries replaced by
        `autosummary` and the module's members added.


    """

    fullfn = os.path.join(app.builder.env.srcdir, fn)

    with io.open(fullfn, encoding='utf8') as fr:
        # Note: we use __name__ here instead of just writing the module name in
        #       case this extension is bundled into another package
        from . import automodapi

        try:
            extensions = app.extensions
        except AttributeError:  # Sphinx <1.6
            extensions = app._extensions

        if automodapi.__name__ in extensions:
            # Must do the automodapi on the source to get the automodsumm
            # that might be in there
            docname = os.path.splitext(fn)[0]
            filestr = automodapi.automodapi_replace(fr.read(), app, True, docname, False)
        else:
            filestr = fr.read()

    spl = _automodsummrex.split(filestr)
    # 0th entry is the stuff before the first automodsumm line
    indent1s = spl[1::5]
    mods = spl[2::5]
    opssecs = spl[3::5]
    indent2s = spl[4::5]
    remainders = spl[5::5]

    # only grab automodsumm sections and convert them to autosummary with the
    # entries for all the public objects
    newlines = []

    # loop over all automodsumms in this document
    for i, (i1, i2, modnm, ops, rem) in enumerate(zip(indent1s, indent2s, mods,
                                                      opssecs, remainders)):
        allindent = i1 + ('    ' if i2 is None else i2)

        # filter out functions-only, classes-only, and ariables-only
        # options if present.
        oplines = ops.split('\n')
        toskip = []
        allowedpkgnms = []
        funcsonly = clssonly = varsonly = False
        for i, ln in reversed(list(enumerate(oplines))):
            if ':functions-only:' in ln:
                funcsonly = True
                del oplines[i]
            if ':classes-only:' in ln:
                clssonly = True
                del oplines[i]
            if ':variables-only:' in ln:
                varsonly = True
                del oplines[i]
            if ':skip:' in ln:
                toskip.extend(_str_list_converter(ln.replace(':skip:', '')))
                del oplines[i]
            if ':allowed-package-names:' in ln:
                allowedpkgnms.extend(_str_list_converter(ln.replace(':allowed-package-names:', '')))
                del oplines[i]
        if [funcsonly, clssonly, varsonly].count(True) > 1:
            msg = ('Defined more than one of functions-only, classes-only, '
                   'and variables-only.  Skipping this directive.')
            lnnum = sum([spl[j].count('\n') for j in range(i * 5 + 1)])
            app.warn('[automodsumm]' + msg, (fn, lnnum))
            continue

        # Use the currentmodule directive so we can just put the local names
        # in the autosummary table.  Note that this doesn't always seem to
        # actually "take" in Sphinx's eyes, so in `Automodsumm.run`, we have to
        # force it internally, as well.
        newlines.extend([i1 + '.. currentmodule:: ' + modnm,
                         '',
                         '.. autosummary::'])
        newlines.extend(oplines)

        ols = True if len(allowedpkgnms) == 0 else allowedpkgnms
        for nm, fqn, obj in zip(*find_mod_objs(modnm, onlylocals=ols)):
            if nm in toskip:
                continue
            if funcsonly and not inspect.isroutine(obj):
                continue
            if clssonly and not inspect.isclass(obj):
                continue
            if varsonly and (inspect.isclass(obj) or inspect.isroutine(obj)):
                continue
            newlines.append(allindent + nm)

    # add one newline at the end of the autosummary block
    newlines.append('')

    return newlines


def generate_automodsumm_docs(lines, srcfn, app=None, suffix='.rst',
                              base_path=None, builder=None,
                              template_dir=None,
                              inherited_members=False):
    """
    This function is adapted from
    `sphinx.ext.autosummary.generate.generate_autosummmary_docs` to
    generate source for the automodsumm directives that should be
    autosummarized. Unlike generate_autosummary_docs, this function is
    called one file at a time.
    """

    from sphinx.jinja2glue import BuiltinTemplateLoader
    from sphinx.ext.autosummary import import_by_name, get_documenter
    from sphinx.util.osutil import ensuredir
    from sphinx.util.inspect import safe_getattr
    from jinja2 import FileSystemLoader, TemplateNotFound
    from jinja2.sandbox import SandboxedEnvironment

    from .utils import find_autosummary_in_lines_for_automodsumm as find_autosummary_in_lines

    if SPHINX_LT_16:
        info = app.info
        warn = app.warn
    else:
        from sphinx.util import logging
        logger = logging.getLogger(__name__)
        info = logger.info
        warn = logger.warning

    # info('[automodsumm] generating automodsumm for: ' + srcfn)

    # Create our own templating environment - here we use Astropy's
    # templates rather than the default autosummary templates, in order to
    # allow docstrings to be shown for methods.
    template_dirs = [os.path.join(os.path.dirname(__file__), 'templates'),
                     os.path.join(base_path, '_templates')]
    if builder is not None:
        # allow the user to override the templates
        template_loader = BuiltinTemplateLoader()
        template_loader.init(builder, dirs=template_dirs)
    else:
        if template_dir:
            template_dirs.insert(0, template_dir)
        template_loader = FileSystemLoader(template_dirs)
    template_env = SandboxedEnvironment(loader=template_loader)

    # read
    # items = find_autosummary_in_files(sources)
    items = find_autosummary_in_lines(lines, filename=srcfn)
    if len(items) > 0:
        msg = '[automodsumm] {1}: found {0} automodsumm entries to generate'
        info(msg.format(len(items), srcfn))

#    gennms = [item[0] for item in items]
#    if len(gennms) > 20:
#        gennms = gennms[:10] + ['...'] + gennms[-10:]
#    info('[automodsumm] generating autosummary for: ' + ', '.join(gennms))

    # remove possible duplicates
    items = list(set(items))

    # keep track of new files
    new_files = []

    # write
    for name, path, template_name, inherited_mem in sorted(items):

        if path is None:
            # The corresponding autosummary:: directive did not have
            # a :toctree: option
            continue

        path = os.path.abspath(os.path.join(base_path, path))
        ensuredir(path)

        try:
            import_by_name_values = import_by_name(name)
        except ImportError as e:
            warn('[automodsumm] failed to import %r: %s' % (name, e))
            continue

        # if block to accommodate Sphinx's v1.2.2 and v1.2.3 respectively
        if len(import_by_name_values) == 3:
            name, obj, parent = import_by_name_values
        elif len(import_by_name_values) == 4:
            name, obj, parent, module_name = import_by_name_values

        fn = os.path.join(path, name + suffix)

        # skip it if it exists
        if os.path.isfile(fn):
            continue

        new_files.append(fn)

        f = open(fn, 'w')

        try:

            if SPHINX_LT_17:
                doc = get_documenter(obj, parent)
            else:
                doc = get_documenter(app, obj, parent)

            if template_name is not None:
                template = template_env.get_template(template_name)
            else:
                tmplstr = 'autosummary_core/%s.rst'
                try:
                    template = template_env.get_template(tmplstr % doc.objtype)
                except TemplateNotFound:
                    template = template_env.get_template(tmplstr % 'base')

            def get_members_mod(obj, typ, include_public=[]):
                """
                typ = None -> all
                """
                items = []
                for name in dir(obj):
                    try:
                        if SPHINX_LT_17:
                            documenter = get_documenter(safe_getattr(obj, name), obj)
                        else:
                            documenter = get_documenter(app, safe_getattr(obj, name), obj)
                    except AttributeError:
                        continue
                    if typ is None or documenter.objtype == typ:
                        items.append(name)
                public = [x for x in items
                          if x in include_public or not x.startswith('_')]
                return public, items

            def get_members_class(obj, typ, include_public=[], exclude_public=[],
                                  include_base=False):
                """
                typ = None -> all
                include_base -> include attrs that are from a base class
                """
                items = []

                # using dir gets all of the attributes, including the elements
                # from the base class, otherwise use __slots__ or __dict__
                if include_base:
                    names = dir(obj)
                else:
                    # Classes deriving from an ABC using the `abc` module will
                    # have an empty `__slots__` attribute in Python 3, unless
                    # other slots were declared along the inheritance chain. If
                    # the ABC-derived class has empty slots, we'll use the
                    # class `__dict__` instead.
                    declares_slots = (
                        hasattr(obj, '__slots__') and
                        not (type(obj) is abc.ABCMeta and
                             len(obj.__slots__) == 0)
                    )

                    if declares_slots:
                        names = tuple(getattr(obj, '__slots__'))
                    else:
                        names = getattr(obj, '__dict__').keys()

                for name in names:
                    try:
                        if SPHINX_LT_17:
                            documenter = get_documenter(safe_getattr(obj, name), obj)
                        else:
                            documenter = get_documenter(app, safe_getattr(obj, name), obj)
                    except AttributeError:
                        continue
                    if typ is None or documenter.objtype == typ:
                        items.append(name)
                public = [x for x in items if x not in exclude_public and
                          (x in include_public or not x.startswith('_'))]
                return public, items

            ns = {}

            if doc.objtype == 'module':
                ns['members'] = get_members_mod(obj, None)
                ns['functions'], ns['all_functions'] = \
                    get_members_mod(obj, 'function')
                ns['classes'], ns['all_classes'] = \
                    get_members_mod(obj, 'class')
                ns['exceptions'], ns['all_exceptions'] = \
                    get_members_mod(obj, 'exception')
            elif doc.objtype == 'class':
                if inherited_mem is not None:
                    # option set in this specifc directive
                    include_base = inherited_mem
                else:
                    # use default value
                    include_base = inherited_members

                ns['members'] = get_members_class(obj, None,
                                                  include_base=include_base)
                ns['methods'], ns['all_methods'] = \
                    get_members_class(obj, 'method', api_class_methods, api_ignore_methods,
                                      include_base=include_base)
                ns['attributes'], ns['all_attributes'] = \
                    get_members_class(obj, 'attribute',
                                      include_base=include_base)
                ns['methods'].sort()
                ns['attributes'].sort()

            parts = name.split('.')
            if doc.objtype in ('method', 'attribute'):
                mod_name = '.'.join(parts[:-2])
                cls_name = parts[-2]
                obj_name = '.'.join(parts[-2:])
                ns['class'] = cls_name
            else:
                mod_name, obj_name = '.'.join(parts[:-1]), parts[-1]

            ns['fullname'] = name
            ns['module'] = mod_name
            ns['objname'] = obj_name
            ns['name'] = parts[-1]

            ns['objtype'] = doc.objtype
            ns['underline'] = len(obj_name) * '='

            # We now check whether a file for reference footnotes exists for
            # the module being documented. We first check if the
            # current module is a file or a directory, as this will give a
            # different path for the reference file. For example, if
            # documenting astropy.wcs then the reference file is at
            # ../wcs/references.txt, while if we are documenting
            # astropy.config.logging_helper (which is at
            # astropy/config/logging_helper.py) then the reference file is set
            # to ../config/references.txt
            if '.' in mod_name:
                mod_name_dir = mod_name.replace('.', '/').split('/', 1)[1]
            else:
                mod_name_dir = mod_name
            if not os.path.isdir(os.path.join(base_path, mod_name_dir)) \
               and os.path.isdir(os.path.join(base_path, mod_name_dir.rsplit('/', 1)[0])):
                mod_name_dir = mod_name_dir.rsplit('/', 1)[0]

            # We then have to check whether it exists, and if so, we pass it
            # to the template.
            if os.path.exists(os.path.join(base_path, mod_name_dir, 'references.txt')):
                # An important subtlety here is that the path we pass in has
                # to be relative to the file being generated, so we have to
                # figure out the right number of '..'s
                ndirsback = path.replace(base_path, '').count('/')
                ref_file_rel_segments = ['..'] * ndirsback
                ref_file_rel_segments.append(mod_name_dir)
                ref_file_rel_segments.append('references.txt')
                ns['referencefile'] = os.path.join(*ref_file_rel_segments)

            rendered = template.render(**ns)
            f.write(cleanup_whitespace(rendered))
        finally:
            f.close()


def setup(app):

    # need autodoc fixes
    # Note: we use __name__ here instead of just writing the module name in
    #       case this extension is bundled into another package
    from . import autodoc_enhancements
    app.setup_extension(autodoc_enhancements.__name__)

    # need inheritance-diagram for automod-diagram
    app.setup_extension('sphinx.ext.inheritance_diagram')

    app.add_directive('automod-diagram', Automoddiagram)
    app.add_directive('automodsumm', Automodsumm)
    app.connect('builder-inited', process_automodsumm_generation)

    app.add_config_value('automodsumm_writereprocessed', False, True)
    app.add_config_value('automodsumm_inherited_members', False, 'env')

    return {'parallel_read_safe': True,
            'parallel_write_safe': True}

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The classes in the astropy docs are documented by their API location,
which is not necessarily where they are defined in the source.  This
causes a problem when certain automated features of the doc build,
such as the inheritance diagrams or the `Bases` list of a class
reference a class by its canonical location rather than its "user"
location.

In the `autodoc-process-docstring` event, a mapping from the actual
name to the API name is maintained.  Later, in the `missing-reference`
event, unresolved references are looked up in this dictionary and
corrected if possible.
"""

from docutils.nodes import literal, reference


def process_docstring(app, what, name, obj, options, lines):
    if isinstance(obj, type):
        env = app.env
        if not hasattr(env, 'class_name_mapping'):
            env.class_name_mapping = {}
        mapping = env.class_name_mapping
        mapping[obj.__module__ + '.' + obj.__name__] = name


def merge_mapping(app, env, docnames, env_other):
    if not hasattr(env_other, 'class_name_mapping'):
        return
    if not hasattr(env, 'class_name_mapping'):
        env.class_name_mapping = {}
    env.class_name_mapping.update(env_other.class_name_mapping)


def missing_reference_handler(app, env, node, contnode):

    if not hasattr(env, 'class_name_mapping'):
        env.class_name_mapping = {}
    mapping = env.class_name_mapping

    reftype = node['reftype']
    reftarget = node['reftarget']
    if reftype in ('obj', 'class', 'exc', 'meth'):
        reftarget = node['reftarget']
        suffix = ''
        if reftarget not in mapping:
            if '.' in reftarget:
                front, suffix = reftarget.rsplit('.', 1)
            else:
                suffix = reftarget

            if suffix.startswith('_') and not suffix.startswith('__'):
                # If this is a reference to a hidden class or method,
                # we can't link to it, but we don't want to have a
                # nitpick warning.
                return node[0].deepcopy()

            if reftype in ('obj', 'meth') and '.' in reftarget:
                if front in mapping:
                    reftarget = front
                    suffix = '.' + suffix

            if (reftype in ('class', ) and '.' in reftarget and
                    reftarget not in mapping):

                if '.' in front:
                    reftarget, _ = front.rsplit('.', 1)
                    suffix = '.' + suffix
                reftarget = reftarget + suffix
                prefix = reftarget.rsplit('.')[0]
                inventory = env.intersphinx_named_inventory
                if (reftarget not in mapping and
                        prefix in inventory):

                    if reftarget in inventory[prefix]['py:class']:
                        newtarget = inventory[prefix]['py:class'][reftarget][2]
                        if not node['refexplicit'] and \
                                '~' not in node.rawsource:
                            contnode = literal(text=reftarget)
                        newnode = reference('', '', internal=True)
                        newnode['reftitle'] = reftarget
                        newnode['refuri'] = newtarget
                        newnode.append(contnode)

                        return newnode

        if reftarget in mapping:
            newtarget = mapping[reftarget] + suffix
            if not node['refexplicit'] and '~' not in node.rawsource:
                contnode = literal(text=newtarget)
            newnode = env.domains['py'].resolve_xref(
                env, node['refdoc'], app.builder, 'class', newtarget,
                node, contnode)
            if newnode is not None:
                newnode['reftitle'] = reftarget
            return newnode


def setup(app):

    app.connect('autodoc-process-docstring', process_docstring)
    app.connect('missing-reference', missing_reference_handler)
    app.connect('env-merge-info', merge_mapping)

    return {'parallel_read_safe': True,
            'parallel_write_safe': True}

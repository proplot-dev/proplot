from docutils import nodes
from os.path import sep
from proplot import rc

# Adapted from matplotlib
# TODO: Understand what the hell is going on here
def rcparam_role(name, rawtext, text, lineno, inliner, options={}, content=[]):
    rendered = nodes.Text(f"rc['{text}']")
    source = inliner.document.attributes['source'].replace(sep, '/')
    relsource = source.split('/docs/', 1)[1]
    levels = relsource.count('/') - 1 # distance to 'docs' folder
    refuri = ('../' * levels + '_build/html/rctools.html?highlight={text}')

    ref = nodes.reference(rawtext, rendered, refuri=refuri)
    node_list = [nodes.literal('', '', ref)]
    if text in rc:
        node_list.append(nodes.Text(f' (default: {rcParamsDefault[text]!r})'))
    return node_list, []


def setup(app):
    app.add_role('rc', rcparam_role)
    return {'parallel_read_safe': True, 'parallel_write_safe': True}

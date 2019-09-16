from docutils import nodes
from os.path import sep
from proplot import rc

# Adapted from matplotlib
# TODO: Understand what the hell is going on here
def rc_role_generator(show_default):
    def rc_role(name, rawtext, text, lineno, inliner, options={}, content=[]):
        rctext = (f"rc['{text}']" if '.' in text else f'rc.{text}')
        rendered = nodes.Text(rctext)
        source = inliner.document.attributes['source'].replace(sep, '/')
        relsource = source.split('/docs/', 1)
        if len(relsource) == 1:
            return [], []
        levels = relsource[1].count('/') # distance to 'docs' folder
        refuri = ('../' * levels + f'en/latest/rctools.html?highlight={text}#'
            + ('rcExtraParams' if '.' in text else 'rcGlobals'))

        ref = nodes.reference(rawtext, rendered, refuri=refuri)
        node_list = [nodes.literal('', '', ref)]
        if show_default and text in rc:
            node_list.append(nodes.Text(' = '))
            node_list.append(nodes.literal('', '', nodes.Text(repr(rc[text]))))
        return node_list, []
    return rc_role

def setup(app):
    app.add_role('rc', rc_role_generator(True))
    app.add_role('rcraw', rc_role_generator(False))
    return {'parallel_read_safe': True, 'parallel_write_safe': True}

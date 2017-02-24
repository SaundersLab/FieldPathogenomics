import nbformat
import re
import os
import time

from fieldpathogenomics.luigi.slurm import SlurmTask


class NotebookTask(SlurmTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def luigi_vars(self, cell_idx):
        '''Operates on cells with a ##luigi-vars magic
           to replace key = value pairs and define new variables'''

        lines = self.nb['cells'][cell_idx]['source'].split('\n')
        for i, l in enumerate(lines):
            if l[0] == '#':
                continue
            try:
                key = re.match('^(\S+)\s*=', l).groups()[0]
                lines[i] = "{0} = {1}".format(key, self.vars_dict[key].__repr__())

            except ValueError:
                raise Exception("Unable to match luigi-vars line {} in {}".format(l, self.notebook))
            except KeyError:
                raise Exception("Unable to match variable {} in {}".format(key, self.notebook))

        self.nb['cells'][cell_idx]['source'] = "\n".join(lines)

    def luigi_meta(self):
        nb_dir, nb_name = os.path.split(self.notebook)
        nb_name = nb_name.split('.ipynb')[0]
        metastring = '\n'.join(["# " + nb_name,
                                "## Created at " + time.strftime("%H:%M:%S on %d/%m/%Y")])

        self.nb['cells'].insert(0, nbformat.v4.new_markdown_cell(metastring))
        self.nb['cells'].insert(0, nbformat.v4.new_code_cell('cd {}'.format(nb_dir)))

    def work(self, *args, **kwargs):

        import nbformat
        import nbconvert
        import re
        import os
        import time

        # Checks out notebook template from repo
        with open(self.notebook, 'r') as fin:
            self.nb = nbformat.read(fin, as_version=4)

        # Populates variable and metadata
        for i, cell in enumerate(self.nb['cells']):
            if re.match('^##luigi-vars', cell['source']):
                self.luigi_vars(i)
        self.luigi_meta()

        # Actually run the notebook here
        ep = nbconvert.preprocessors.ExecutePreprocessor()
        ep.preprocess(self.nb, {})

        # Make the HTML conversion
        html = nbconvert.html.HTMLExporter()
        (body, resources) = html.from_notebook_node(self.nb)

        # Make use of luigi's built in atomic file writing
        with open(self.output().path.split('.ipynb')[0] + '.html', 'w') as fout:
            fout.write(body)

        with self.output().open('w') as fout:
            nbformat.write(self.nb, fout)

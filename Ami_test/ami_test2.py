import utils_func as uf
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--conf", default = "user.osalin.MadGraph_WmZ_llqq_FM1_QUAD")
opts, _ = parser.parse_args()

EFT_op, proces, decay= uf.extract_EFT_op_proces_dec(opts.conf)
xsection_fb = uf.cross_section_fb(EFT_op, proces, decay)
print(f'Cross section in fb: {xsection_fb}')


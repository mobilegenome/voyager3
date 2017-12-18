from dna_features_viewer import BiopythonTranslator
from os import listdir
from os.path import join

filepath = "/home/fritjof/fin_whale_genome_assembly/data/2017-11-01_L1seq/annotation/LATEST/"
gb_files = [join(filepath, f) for f in listdir(filepath)]
print(gb_files)

graphic_record = BiopythonTranslator().translate_record(gb_files[0])

ax,_ = graphic_record.plot(figure_width=20)
ax.figure
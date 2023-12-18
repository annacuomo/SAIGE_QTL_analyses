import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# example
venn2(subsets = (10, 5, 2), set_labels = ('Group A', 'Group B'))
plt.show()

# below is super manual, using numbers from genes_overlap.R script
venn2(subsets = (887, 151, 3103), set_labels = ('SAIGE-QTL', 'TensorQTL'), set_colors=('#882E72','pink'))
plt.title("CD4 NC")
plt.show()

venn2(subsets = (312, 50, 793), set_labels = ('SAIGE-QTL', 'TensorQTL'), set_colors=('#B178A6','pink'))
plt.title("CD4 ET")
plt.show()

venn2(subsets = (19, 5, 34), set_labels = ('SAIGE-QTL', 'TensorQTL'), set_colors=('#D6C1DE','pink'))
plt.title("CD4 SOX4")
plt.show()

venn2(subsets = (975, 72, 1525), set_labels = ('SAIGE-QTL', 'TensorQTL'), set_colors=('#1965B0','pink'))
plt.title("CD8 ET")
plt.show()

venn2(subsets = (524, 87, 1363), set_labels = ('SAIGE-QTL', 'TensorQTL'), set_colors=('#5289C7','pink'))
plt.title("CD8 NC")
plt.show()

venn2(subsets = (308, 26, 346), set_labels = ('SAIGE-QTL', 'TensorQTL'), set_colors=('#7BAFDE','pink'))
plt.title("CD8 S100B")
plt.show()

venn2(subsets = (682, 111, 1658), set_labels = ('SAIGE-QTL', 'TensorQTL'), set_colors=('#4EB265','pink'))
plt.title("NK")
plt.show()

venn2(subsets = (93, 11, 133), set_labels = ('SAIGE-QTL', 'TensorQTL'), set_colors=('#90C987','pink'))
plt.title("NK R")
plt.show()

venn2(subsets = (56, 4, 57), set_labels = ('SAIGE-QTL', 'TensorQTL'), set_colors=('#CAE0AB','pink'))
plt.title("Plasma")
plt.show()

venn2(subsets = (367, 46, 550), set_labels = ('SAIGE-QTL', 'TensorQTL'), set_colors=('#F7EE55','pink'))
plt.title("B Mem")
plt.show()

venn2(subsets = (440, 36, 681), set_labels = ('SAIGE-QTL', 'TensorQTL'), set_colors=('#F6C141','pink'))
plt.title("B IN")
plt.show()

venn2(subsets = (661, 14, 276), set_labels = ('SAIGE-QTL', 'TensorQTL'), set_colors=('#F1932D','pink'))
plt.title("Mono C")
plt.show()

venn2(subsets = (399, 18, 269), set_labels = ('SAIGE-QTL', 'TensorQTL'), set_colors=('#E8601C','pink'))
plt.title("Mono NC")
plt.show()

venn2(subsets = (119, 17, 149), set_labels = ('SAIGE-QTL', 'TensorQTL'), set_colors=('#DC050C','pink'))
plt.title("DC")
plt.show()

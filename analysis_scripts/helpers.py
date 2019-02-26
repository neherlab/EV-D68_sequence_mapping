def name_translations(fname):
	name_translation_table = {}
	with open(fname) as fh:
		for line in fh:
			labid, strain = line.strip().split('\t')
			name_translation_table[labid] = strain
	return name_translation_table


def add_panel_label(ax, label, fs, x_offset=-0.1):
    '''Add a label letter to a panel'''
    ax.text(x_offset, 0.95, label,
            transform=ax.transAxes,
            fontsize=fs*1.5)



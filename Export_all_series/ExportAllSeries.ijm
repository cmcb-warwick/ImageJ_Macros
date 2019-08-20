
run("Bio-Formats Macro Extensions");
id = File.openDialog("Choose a file");
Ext.setId(id);
Ext.getSeriesCount(seriesCount);
for (s=0; s<seriesCount; s++) {
  
  run("Bio-Formats Importer",
		"open=[" + id + "] " +
		"autoscale " +
		"color_mode=Grayscale " +
		"view=Hyperstack " +
		"stack_order=XYCZT " +
		"series_" + (s+1));
  run("Bio-Formats Exporter", "save=["+File.getParent(id)+"/"+File.nameWithoutExtension+"_series_"+s+".ome.tiff] compression=Uncompressed");
  close();
}
Ext.close();
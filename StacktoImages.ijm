/*
 * Macro template to process multiple images in a folder
 */

#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".tif") suffix
setBatchMode(true); 

print(input);
print(output);
processFolder(input);

function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	print(list.length);
	
	
	for (imcnt = 0; imcnt < list.length; imcnt++) {
	//for (imcnt = 0; imcnt < 6; imcnt++) {
		if(File.isDirectory(input + File.separator + list[imcnt]))
			processFolder(input + File.separator + list[imcnt]);
			print(list[imcnt]);
		//if(imcnt>0){
		if(endsWith(list[imcnt], suffix))
			processFile(input, output, list[imcnt]);
		//}
	}

function processFile(input, output, file) {

	print("Processing: " + input + File.separator + file);
	print("Saving to: " + output);
	close("*");

	path = input + File.separator + file;
	dir = File.getParent(path);
	name = File.getName(path);

	
	open(path);
	run("Stack to Images");
	outputPath = output + File.separator + file;
	run("Image Sequence... ", "format=TIFF save=["+outputPath+"]");
	close();
	run("Image Sequence... ", "format=TIFF save=["+outputPath+"]");
	close();
	run("Image Sequence... ", "format=TIFF save=["+outputPath+"]");
	close();
}
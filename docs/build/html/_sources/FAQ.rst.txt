FAQ
===

* Can I run ConsensusCruncher in parallel?
	Yes, we provide ``generate_scripts.sh`` as a wrapper to create scripts for each 
	individual sample. If you have access to a cluster, you can submit these scripts 
	as separate jobs to be run in parallel. 

* Why is ConsensusCruncher taking so long? 
	The time it takes ConsensusCruncher to run depends on the sequencing depth, number of 
	UMIs, and the computer you're using. Time can range from a few minutes to several hours 
	(& up to a day with some of the larger samples). Please be patient!
	
* Why is ConsensusCruncher taking so much memory?
	The ``consensus`` mode of ConsensusCruncher groups reads of a BAM file by UMI and 
	condenses it into consensus sequences. The amount of memory is dependent on the size of
	your BAM files. If you have large BAMs and have access to a cluster, we recommend 
	running ConsensusCruncher on the high memory queue.

* Why do we need a cytoband bedfile to separate the data?
	Depending on the size of your BAM file, it may be too large to load all the reads into
	memory. We use coordinates of cytobands to split our data, but you could also use a 
	target bedfile if you're only interested in specific regions of the genome. 
	
* Can I input BAMs into ConsensusCruncher?
	Yes, as long as your BAM files have UMIs in the queryname of each read, it can be used
	for the ``consensus`` mode of ConsensusCruncher. (Please note: UMIs need to be separated
	by '|' and paired UMIs split by '.'). e.g. "HWI-D00331:196:C900FANXX:5:1101:11551:2948|AC.CT"

* Can I use ConsensusCruncher for single-end sequencing data?
	Yes, you can use the ``SSCS_maker.py`` directly to generate consensus sequences. However,
	please note that single-end sequencing data won't benefit from Singleton Correction
	as that requires paired-end data. 
	
Questions
---------
If you have any questions, please post it as an issue on our `GitHub <https://github.com/pughlab/ConsensusCruncher/blob/master/README.md>`_.

Alternatively, you can contact:
Nina Wang (nina.tt.wang@gmail.com), Trevor Pugh (Trevor.Pugh@uhn.ca), Scott Bratman (Scott.Bratman@rmp.uhn.ca)

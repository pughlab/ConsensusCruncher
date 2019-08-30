How it works
============

.. image:: https://user-images.githubusercontent.com/13406244/63974147-468d5880-ca7a-11e9-849d-289e69601771.png
   :width: 600

Unique molecular identifiers (UMIs) composed of molecular barcodes and sequence features are used aggregate reads derived from the same strand of a template molecule. Amalgamation of such reads into single strand consensus sequences (SSCS) removes discordant bases, which effectively eliminates polymerase and sequencer errors. Complementary SSCSs can be subsequently combined to form a duplex consensus sequence (DCS), which eliminates asymmetric strand artefacts such as those that develop from oxidative damage.

Conventional UMI-based strategies rely on redundant sequencing from both template strands to form consensus sequences and cannot error suppress single reads (singleton). We enable singleton correction using complementary duplex reads in the absence of redundant sequencing.

**ConsensusCruncher schematic:**

- An uncollapsed bamfile is first processed through SSCS_maker.py to create an 
  error-suppressed single-strand consensus sequence (SSCS) bamfile and an uncorrected 
  singleton bamfile.	
- The singletons can be corrected through singleton_correction.py, which error suppress 
  singletons with its complementary SSCS or singleton read.
- SSCS reads can be directly made into duplex consensus sequences (DCS) or merged with 
  corrected singletons to create an expanded pool of DCS reads (Figure illustrates singleton 
  correction merged work flow).
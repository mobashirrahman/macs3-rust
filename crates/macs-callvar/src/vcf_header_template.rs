/// The VCF header template, transcribed verbatim from `callvar_cmd.py`'s
/// `VCFHEADER`.
///
/// `callvar_cmd.py` defines **two** templates -- `VCFHEADER_0` and `VCFHEADER` --
/// and `run()` uses `VCFHEADER` (line 162). They differ: `VCFHEADER` declares the
/// `DBIC` INFO field that `VCFHEADER_0` lacks, and orders the genotype-likelihood
/// fields differently. Reading the wrong one is the easy mistake here, and it shifts
/// every `##` line from 13 onwards.
///
/// `%s` substitutions are applied by [`vcf_header`]:
///   0. the run date as `%Y%m%d`,
///   1. the MACS version string,
///   2. the `Program_Args` value (see [`program_args`]).
pub const VCF_HEADER_TEMPLATE: &str = r##"##fileformat=VCFv4.1
##fileDate=%s
##source=MACS_V%s
##Program_Args=%s
##INFO=<ID=M,Number=.,Type=String,Description="MACS Model with minimum BIC value">
##INFO=<ID=MT,Number=.,Type=String,Description="Mutation type: SNV/Insertion/Deletion">
##INFO=<ID=DPT,Number=1,Type=Integer,Description="Depth Treatment: Read depth in ChIP-seq data">
##INFO=<ID=DPC,Number=1,Type=Integer,Description="Depth Control: Read depth in control data">
##INFO=<ID=DP1T,Number=.,Type=String,Description="Read depth of top1 allele in ChIP-seq data">
##INFO=<ID=DP2T,Number=.,Type=String,Description="Read depth of top2 allele in ChIP-seq data">
##INFO=<ID=DP1C,Number=.,Type=String,Description="Read depth of top1 allele in control data">
##INFO=<ID=DP2C,Number=.,Type=String,Description="Read depth of top2 allele in control data">
##INFO=<ID=DBIC,Number=.,Type=Float,Description="Difference of BIC of selected model vs second best alternative model">
##INFO=<ID=BICHOMOMAJOR,Number=1,Type=Integer,Description="BIC of homozygous with major allele model">
##INFO=<ID=BICHOMOMINOR,Number=1,Type=Integer,Description="BIC of homozygous with minor allele model">
##INFO=<ID=BICHETERNOAS,Number=1,Type=Integer,Description="BIC of heterozygous with no allele-specific model">
##INFO=<ID=BICHETERAS,Number=1,Type=Integer,Description="BIC of heterozygous with allele-specific model">
##INFO=<ID=AR,Number=1,Type=Float,Description="Estimated allele ratio of heterozygous with allele-specific model">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth after filtering bad reads">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality score">
##FORMAT=<ID=PL,Number=3,Type=Integer,Description="Normalized, Phred-scaled genotype likelihoods for 00, 01, 11 genotype">"##;

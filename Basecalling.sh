#Install dorado throug conda, then run the following

dorado basecaller hac path/to/pod5 > output.bam

#Next, we can trim the sequencing adapters

dorado trim output.bam --sequencing_kit sequencing_kit > trimmed.bam
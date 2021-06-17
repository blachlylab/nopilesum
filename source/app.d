import std.stdio;
import std.getopt;
import std.algorithm : filter, count, map;
import std.range : repeat;
import std.array : array;
import std.format : format;

import dhtslib.vcf;
import dhtslib.sam;
import dhtslib.coordinates;
import htslib.hts_log;

import filters;

float maxAF = 0.2;
float minAF = 0.01;
bool verbose = false;
bool verbose2 = false;

int main(string[] args)
{
	// parse some args
	auto res = getopt(args, 
			"max-af","VCF records with INFO/AF fields above this threshold won't be used (default 0.2)", &maxAF, 
			"min-af", "VCF records with INFO/AF fields below this threshold won't be used (default 0.01)", &minAF,
			"verbose|v", "see warnings", &verbose,
			"debug", "see extra info", &verbose2
		);
	
	// set logging
	hts_set_log_level(htsLogLevel.HTS_LOG_ERROR);
	if(verbose) hts_set_log_level(htsLogLevel.HTS_LOG_WARNING);
	if(verbose2) hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
	if (res.helpWanted || (args.length < 3))
	{
		defaultGetoptPrinter("nopilesum: get pileup summaries\nusage: nopilesum <in.bam> <in.vcf> > summary.txt",
				res.options);
		stderr.writeln();
		return 0;
	}
	
	// set up readers
	auto bamr = SAMReader(args[1]);
	auto vcfr = VCFReader(args[2]);

	// get sample info if exists
	if(RecordType.RG in bamr.header){
		auto sample = bamr.header.valueByPos(RecordType.RG, 0, "SM").idup;
		writeln("#<METADATA>SAMPLE=" ~ sample);
	}
	
	// write header
	writefln("%s\t%s\t%s\t%s\t%s\t%s","contig","position","ref_count","alt_count","other_alt_count","allele_frequency");
	foreach (VCFRecord rec; vcfr)
	{
		// get AF field from vcf
        auto AF = getGermlineResourceAF(rec);
		// make sure AF is within range
        if(AF > maxAF || AF < minAF){
            hts_log_info("nopilesum:filters", "VCFRecord AF INFO is either greater than max or less than min");
            continue;
        }
		// query bam for this pos and count alleles
		auto alleles =  bamr.query(rec.chrom, rec.coordinates).countAlleles(rec.pos);

		// just some debug code
        //if(bamr.query(rec.chrom, rec.coordinates).count > 0){
        //    writefln("%s\t%d\t%d\t%d",
        //            rec.chrom,
        //            rec.pos + 1,
        //            bamr.query(rec.chrom, rec.coordinates).count,
        //            bamr.query(rec.chrom, rec.coordinates).filter!passesFilters.filter!(x => rec.pos >= x.pos).count);
        //    writeln(alleles," ", rec.altAllelesAsArray[0]);
        //}

		int refCount;
		int altCount;
		int otherAltCount;
		string alt = rec.altAllelesAsArray[0];

		// if a complex INDEL we ignore, to difficult to handle
		// Does GATK handle this?
		if(rec.refAllele.length > 1 && alt.length > 1){
			hts_log_warning("nopilesum:main", "we do not handle complex indels, skipping...");
			continue;
		}else if(rec.refAllele.length > 1 && alt.length == 1){
			alt = '-'.repeat(rec.refAllele.length).array.idup;
		}
		// loop over alleles and tally
		// dels are a special case they are represented as dashes of
		// the same length of the deletion
		// this keeps us from needing an MD string for getAligned pairs
		foreach (allele; alleles.byKey)
		{
			if(rec.refAllele == allele){
				refCount = alleles[rec.refAllele];
			}
			else if(alt == allele){
				altCount = alleles[alt];
			}else if(allele == "N"){
                continue;   
            }else{
				otherAltCount += alleles[allele];
			}

		}
		// if any depth report 
		if(refCount + altCount + otherAltCount != 0){
            //writeln(alleles);
			writefln("%s\t%d\t%d\t%d\t%d\t%s",
						rec.chrom, 
						rec.pos + 1,
                        refCount,
						altCount,
						otherAltCount,
						AF);
		}
	}
	return 0;
}

/// loop over all reads and tally alleles at the singular position
/// dels are a special case they are represented as dashes of
/// the same length of the deletion
/// this keeps us from needing an MD string for getAligned pairs
int[string] countAlleles(SAMRange)(SAMRange records, ZeroBased pos)
{
	int[string] alleles;
	foreach (SAMRecord rec; records.filter!passesFilters.filter!(x => pos >= x.pos))
	{
		auto pairs = rec.getAlignedPairs!false(pos - rec.pos, (pos - rec.pos) + 1);

		// if SNP, just need the first aligned pair
		if(pairs.front.cigar_op == Ops.MATCH || pairs.front.cigar_op == Ops.EQUAL || pairs.front.cigar_op == Ops.DIFF)
		{
			alleles[[pairs.front.queryBase].idup]++;
		
		// if a DEL, add an allele of dashes the same length as the deletion
		}else if(pairs.front.cigar_op == Ops.DEL){
			alleles['-'.repeat(pairs.count).array.idup]++;
		// if a INS, add the whole insertion
		}else if(pairs.front.cigar_op == Ops.INS){
			alleles[pairs.map!(x => x.queryBase).array.idup]++;
		// extra case
		}else{
			hts_log_warning("nopilesum:countAlleles", "saw strange cigar op: %s in read %s".format(CIGAR_STR[pairs.front.cigar_op], rec.queryName));
			continue;
		}
	}
	return alleles;
}

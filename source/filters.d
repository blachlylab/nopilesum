module filters;

import std.algorithm : canFind, sum, map;
import std.format : format;
import std.utf : toUTFz;

import dhtslib.sam.record;
import dhtslib.sam.cigar;
import dhtslib.vcf.record;
import dhtslib.vcf.header;
import htslib.vcf;
import htslib.hts_log;
import htslib.sam;

pragma(inline, true):
/// does alignment pass all filters
bool passesFilters(SAMRecord rec)
{
    return rec.isQCPass && 
        rec.isMapped &&
        rec.isNotDuplicate && 
        rec.isPrimary && 
        rec.hasGoodMate && 
        rec.hasPositiveQual &&
        rec.hasGoodCigar &&
        rec.hasPositiveAlignedLength;
}


/// is cigar generally well formed?
bool hasGoodCigar(SAMRecord rec){
    Cigar cigar = rec.cigar;
    bool ret = true;
    /// check for starting del
    if(cigar[0].isClipping)
        cigar = cigar[1..$];
    if(cigar[0].isClipping)
        cigar = cigar[1..$];
    if(cigar[0].op == Ops.DEL)
        ret = false;

    /// check for ending del
    if(cigar[$-1].isClipping)
        cigar = cigar[0..$-1];
    if(cigar[$-1].isClipping)
        cigar = cigar[0..$-1];
    if(cigar[$-1].op == Ops.DEL)
        ret = false;

    /// check for mapped?
    if(cigar.alignedLength <= 0) ret = false;

    /// check for consecutive INDEL ops
    Ops last_op;
    foreach (Ops op; CigarItr(cigar))
    {
        if(last_op == Ops.DEL && op == Ops.INS)
            ret = false;
        if(last_op == Ops.INS && op == Ops.DEL)
            ret = false;
        last_op = op;
    }
    debug if (!ret) hts_log_warning("nopilesum:filters", "hasGoodCigar filter tripped by read %s".format(rec.queryName));
    return ret;
}

/// does cigar have positive aligned length?
bool hasPositiveAlignedLength(SAMRecord rec)
{
    /// check for invalid cigar
    bool ret = rec.cigar.alignedLength > 0 ? true : false;
    debug if (!ret) hts_log_warning("nopilesum:filters", "hasPositiveAlignedLength filter tripped by read %s".format(rec.queryName));
    return ret;
}

/// does read pass QC?
bool isQCPass(SAMRecord rec)
{
    bool ret = !(cast(bool)(rec.b.core.flag & BAM_FQCFAIL));
    debug if (!ret) hts_log_warning("nopilesum:filters", "isQCPass filter tripped by read %s".format(rec.queryName));
    return ret;
}

/// is read a Not a duplicate
bool isNotDuplicate(SAMRecord rec)
{
    bool ret = !(cast(bool)(rec.b.core.flag & BAM_FDUP));
    debug if (!ret) hts_log_warning("nopilesum:filters", "isNotDuplicate filter tripped by read %s".format(rec.queryName));
    return ret;
}

/// is primary alignment?
bool isPrimary(SAMRecord rec)
{
    bool ret = !rec.isSecondary && !rec.isSupplementary;
    debug if (!ret) hts_log_warning("nopilesum:filters", "isPrimary filter tripped by read %s".format(rec.queryName));
    return ret;
}

/// is mate good?
bool hasGoodMate(SAMRecord rec)
{
    bool ret = !rec.isMateMapped || rec.tid == rec.mateTID;
    debug if (!ret) hts_log_warning("nopilesum:filters", "hasGoodMate filter tripped by read %s".format(rec.queryName));
    return ret;
}

/// Mapping qual is not zero 
bool hasPositiveQual(SAMRecord rec)
{
    bool ret = rec.qual > 0;
    debug if (!ret) hts_log_warning("nopilesum:filters", "hasPositiveQual filter tripped by read %s".format(rec.queryName));
    return ret;
}

/// read is well formed
bool isReadWellFormed(SAMRecord rec)
{
    
    auto ret = rec.pos > 0 && 
        ((rec.pos + rec.cigar.alignedLength) > rec.pos) &&
        !(CigarItr(rec.cigar).canFind!(x => x == Ops.REF_SKIP)) && 
        rec.length > 0 &&
        !(CigarItr(rec.cigar).map!(x => cast(int) (CigarOp(0,x).is_query_consuming)).sum == rec.length);
    debug if (!ret) hts_log_warning("nopilesum:filters", "isReadWellFormed filter tripped by read %s".format(rec.queryName));
    return ret; 
}

/// return true if VCF record AF field is <= max and >= min 
bool filterGermlineResourceVCFRecords(VCFRecord rec, float max, float min)
{
    bcf1_t * b = rec.line;
    bcf_hdr_t * h = rec.vcfheader.hdr;
    auto afRec = bcf_get_info(h, b, toUTFz!(char *)("AF"));
    if(!afRec){
        debug hts_log_warning("nopilesum:filters", "VCFRecord doesn't contian AF INFO field");
        return false;
    }
    if(afRec.type != BCF_BT_FLOAT){
        debug hts_log_warning("nopilesum:filters", "VCFRecord AF INFO field isn't a float");
        return false;
    }
    if(afRec.v1.f > max || afRec.v1.f < min){
        debug hts_log_debug("nopilesum:filters", "VCFRecord AF INFO is either greater than max or less than min");
        return false;
    }
    return true;
}

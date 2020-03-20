#!/usr/bin/python3
import sys
import timeit

class Exon:

    def __init__(self, tr, symbol, num, c, start, end, strand, islast=True):
        self.id = tr+"_"+str(num)
        self.transcript = tr
        self.symbol = symbol
        self.num = num
        self.chr = c
        self.start = start
        self.end = end
        self.strand = strand
        self.length = end - start + 1 # start & end are 1-based
        self.depths = []
        self.islast = islast

    def add_depth(self, d):
        self.depths.append(d)

    def total_coverage(self):
        return sum(self.depths)

    def pct_above(self, t):
        return (len([d for d in self.depths if d > t])*100) / self.length

    def __str__(self):
        return "\t".join([self.chr, str(self.start), str(self.end), self.strand, self.transcript, str(self.num), self.symbol])


class Transcript:
    def __init__(self, id, gene, exons):
        self.id = id
        self.gene = gene
        self.exons = exons
        self.chr = exons[0].chr

    @property
    def start(self):
        return exons[0].start

    @property
    def end(self):
        return exons[-1].end

    @property
    def length(self):
        return sum([e.length for e in self.exons])

    def total_coverage(self):
        return sum([e.total_coverage() for e in self.exons])

    def pct_above(self, t):
        cnt = 0
        for e in self.exons:
            cnt += len([d for d in e.depths if d >= t])
        return (cnt*100)/self.length

    def append_exon(self, exon):
        self.exons[-1].islast = False
        self.exons.append(exon)


def add2all(exons, depth):
    for e in exons:
        e.add_depth(depth)


def write_exon_stats(exon, exons_stats_stream):
    rec = '\t'.join([exon.id, exon.transcript, exon.symbol]+[str(round(e,2)) for e in [exon.total_coverage()/exon.length, exon.pct_above(10), exon.pct_above(20), exon.pct_above(30)]])
    exons_stats_stream.write(rec+'\n')


def write_transcript_stats(tr, transcript_stats_stream):
    rec = '\t'.join([tr.id, tr.gene]+[str(round(e,2)) for e in [tr.total_coverage()/tr.length, tr.pct_above(10), tr.pct_above(20), tr.pct_above(30)]])
    transcript_stats_stream.write(rec + '\n')


#def get_closest_end(exons):
#    closest = exons[0].end
#    for i,e in enumerate(exons):
#        if e.start > closest:
#            return closest
#        if e.end < closest:
#            closest = e.end

def insert_ordered_by_end(exons, exon):
    for i in range(len(exons),0,-1):
        if exons[i-1].end <= exon.end:
            return exons[0:i] + [exon] + exons[i:]
    return [exon]+exons

def write_stats(depth_stream, exons, transcripts, exons_stats_stream, transcript_stats_stream):

    # exons and transcripts with start<=curr_pos and end>=curr_pos, sorted by end
    max_curr_exons=0
    block_size=10000000
    block_time_start=0

    curr_exons = []

    last_contig=''
    next_start_c = exons[0].chr
    next_start_p = exons[0].start
    next_end_c = '' # this is updated when exon is added to curr_exons
    next_end_p = 0  # this is updated when exon is added to curr_exons

#    print(f"Next start {next_start_c}:{next_start_p}")
#    print(f"Next end {next_end_c}:{next_end_p}")

    for l in depth_stream:
        ls = l.strip().split('\t')
        contig, pos, depth = ls[0], int(ls[1]), int(ls[2])

        if contig != last_contig:
            if last_contig != '':
                print(f"###################\n"+\
                      f"End of contig\n"+\
                      f"Current exons: {len(curr_exons)}\n"+\
                      f"Max current exons: {max_curr_exons}\n"+\
                      f"Time for block (partial): {round(timeit.default_timer() - block_time_start,2)}s.\n")
            print(f"Starting contig {contig}")
            last_contig = contig
            block_time_start = timeit.default_timer()
            max_curr_exons=0
           
        if contig==next_start_c and pos == next_start_p:
            # reached new exon - add it to current
#            print(f"Position {contig}:{pos}")
            while len(exons)>0 and exons[0].chr == contig and exons[0].start == pos:
                exon = exons.pop(0)
#                print(f"New exon: {exon}")
                curr_exons = insert_ordered_by_end(curr_exons, exon)

            if len(exons)>0:
                next_start_c,next_start_p = exons[0].chr,exons[0].start
                next_end_c,next_end_p = curr_exons[0].chr,curr_exons[0].end
            else:  # do nothing as there is no features to work on and the next_start will not happen again
                pass
            
#            print(f"Current exons: {len(curr_exons)}. Next start {next_start_c}:{next_start_p}")


        add2all(curr_exons, depth)

        if contig == next_end_c and pos == next_end_p:
            # reached end of an exon - remove it from current and write it out
#            print(f"Position {contig}:{pos}")
            while len(curr_exons)>0 and curr_exons[0].chr == contig and curr_exons[0].end == pos:
                exon = curr_exons.pop(0)
#                print(f"Removing exon: {exon}")
                write_exon_stats(exon, exons_stats_stream)
                if exon.islast:    # last exon of transcript
#                    print(f"Removing transcript: {exon.transcript}")
                    write_transcript_stats(transcripts[exon.transcript], transcript_stats_stream)
                    del transcripts[exon.transcript]

            if len(curr_exons)>0:  # if we are still within an exon
                next_end_c,next_end_p = curr_exons[0].chr,curr_exons[0].end
            else:  # we are outside of any exons, positions will get updated when next exon is added into curr_exons
                pass #next_end_c,next_end_p = exons[0].chr,get_closest_end(exons)

#            print(f"Current exons: {len(curr_exons)}. Next end {next_end_c}:{next_end_p}")

        max_curr_exons = max(max_curr_exons, len(curr_exons))
        if pos % block_size == 0:
            print(f"###################\n"+\
                  f"Position: {contig}:{pos}\n"+\
                  f"Current exons: {len(curr_exons)}\n"+\
                  f"Max current exons: {max_curr_exons}\n"+\
                  f"Time per {block_size/1000000}Mb: {round(timeit.default_timer() - block_time_start,2)}s.\n")
            block_time_start = timeit.default_timer()
            max_curr_exons=0


#def insert_ordered_by_start(exon_list, exons_to_insert):
#    for exon in exons_to_insert:
#        i = -1
#        for i in range(len(exon_list),0,-1):
#            if exon_list[i-1].start <= exon.start:
#                break
#        exon_list = exon_list[0:i]+[exon]+exon_list[i:] if i > 0 else [exon]+exon_list
#
#    return exon_list



def read_exon_map(fname):
    cols = {"CHR": 0, "START": 1, "END": 2, "STRAND": 3, 
            "TR_ID": 4, "EXON_NUM": 5, "GENE_SYMBOL": 6}
                
    exons=[]
    transcripts={}
    with open(fname) as f:
        for line in f:
            ls = line.strip().split("\t")
            tr_id = ls[cols["TR_ID"]]
            gene_symbol = ls[cols["GENE_SYMBOL"]]
            # add 1 to start because BED and UCSC tables are 0-based, while samtools depth output is 1-based
            exon = Exon(tr_id, gene_symbol, int(ls[cols["EXON_NUM"]]), 
                        ls[cols["CHR"]], int(ls[cols["START"]])+1, int(ls[cols["END"]]), ls[cols["STRAND"]])
            exons.append(exon)
            if tr_id in transcripts:
                transcripts[tr_id].append_exon(exon)
            else:
                transcripts[tr_id] = Transcript(tr_id, gene_symbol, [exon])

    return exons, transcripts


def make_exons(tr_id, symbol, contig, starts, ends, strand):
    s = starts.strip(',').split(',')
    e = ends.strip(',').split(',')
    exons = []
    for i in range(0, len(s)):
        exons.append(Exon(tr_id, symbol, i+1, contig, int(s[i]), int(e[i]), strand))
    return exons


def split_transcript_map_into_exon_map(transcript_stream, exon_stream):
    # bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames
    cols = {"TR_ID": 1, "GENE_SYMBOL": 12,
                "CHR": 2, "STRAND":3, "TR_START": 4, "TR_END": 5,
                "CDS_START": 6, "CDS_END": 7,
                "EXON_NUM": 8, "EXON_STARTS": 9, "EXON_ENDS": 10}

    for line in transcript_stream:
        if line.startswith("#"):
            continue
        ls = line.strip().split("\t")
        es = make_exons(ls[cols["TR_ID"]], ls[cols["GENE_SYMBOL"]], 
                        ls[cols["CHR"]], ls[cols["EXON_STARTS"]], ls[cols["EXON_ENDS"]], ls[cols["STRAND"]])
        for e in es:
            exon_stream.write(str(e)+'\n')



if __name__ == "__main__":
    """
	Usage:
	1) Download NCBIrefSeq transcript table
		refseq.txt
        2) Remove non-standard chromosomes and sort by chr and position, e.g.
		grep -vP 'chr[0-9a-zA-Z]+_' refseq.txt | sort -k3,3V -k5,5n -k6,6n > refseq.sorted.txt
	3) Extract exons and sort by coordinate
		cat refseq.sorted.txt | ./gene_depth.py | bedtools sort -i - > refseq.exons.bed
        4) Calculate depth
		samtools depth -a BAM | ./gene_depth.py refseq.exons.bed OUTPUT_PREFIX
	5) Mean, median, and fraction of bp covered by 10, 20 and 30x will be provided in two files:
		OUTPUT_PREFIX.exon_stats
		OUTPUT_PREFIX.transcript_stats
	
	TODOs
	Stats for transcripts present on two chromosomes (e.g. X and Y) are not separated and may be wrong

    """


    if len(sys.argv)==1:  
        # conversion of the refseq transcript table into exon table
        # input on stdin, output on stdout
        
        split_transcript_map_into_exon_map(sys.stdin, sys.stdout)

    else:
        # depth calculation
        # samtools depth output on stdin, exon table and prefix for output as arguments

        depth_stream = sys.stdin

        print("Loading exon coordinates...")
        exons,transcripts = read_exon_map(sys.argv[1])

        print("Done. Calculating coverage...")
        with open(sys.argv[2]+"_exonstats", 'w') as ex_stats, \
             open(sys.argv[2]+"_transcriptstats", 'w') as tr_stats:
             write_stats(depth_stream, exons, transcripts,
                         ex_stats, tr_stats)



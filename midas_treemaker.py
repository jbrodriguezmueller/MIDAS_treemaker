#!/usr/bin/python
from Bio import SeqIO, AlignIO, Phylo
import os, argparse
import subprocess
from collections import defaultdict

def read_many_alns( list_of_files, debug=False ):
  """Reads a list of (fasta aln) filenames, returns dict with headers:seqs"""
  all_seqs = {}
  for filex in list_of_files:
    fasta_sequences = AlignIO.parse(open(filex),'fasta')
    fasta_list = fasta_sequences.next()
    if debug:
      print "@ read_many_alns, fasta_list is "
      print fasta_list
    for fasta in fasta_list:
      name, sequence = fasta.id, fasta.seq
      all_seqs[ name ] = sequence
  if debug:
    print("Example sequence: ")
    print(sequence)
  return all_seqs

def nominal_lens( all_seqs ):
  all_lens = defaultdict( list )
  for key in all_seqs:
    pieces = key.split("_")
    all_lens[ pieces[-1] ].append( len( all_seqs[key] ) )
  for key in all_lens:
    print("Working on "+key)
    print("It has "+str(len(all_lens[key]))+" elements, ")
    #print all_seqs[all_seqs.keys()[0]]
    newstuff = len( set( all_lens[key]))
    if newstuff > 1:
      print "ERROR!"
      print "all_lens = "
      print all_lens
      exit()
    else:
      all_lens[key] = all_lens[key][0]
  return all_lens

def read_ids( input_file ):
  """Reads a file with ids, one per line, # at position 0 => comment"""
  in_file = open( input_file )
  in_data = []
  for rawline in in_file:
    if rawline[0] != "#":
      in_data.append( rawline.strip() )
  in_file.close()
  return in_data

def query_species( all_seqs, species ):
  sp_matches = {}
  for key in all_seqs:
    pieces = key.split("_")
    if species in pieces[0]:
      sp_matches[ key ] = all_seqs[key]
  return sp_matches
def query_protein( all_seqs, protein ):
  sp_matches = {}
  for key in all_seqs:
    pieces = key.split("_")
    if protein in pieces[-1]:
      sp_matches[ key ] = all_seqs[key]
  return sp_matches

def get_species_from_db( argsdict ):
  fasta_sequences = SeqIO.parse(open(argsdict['df']),'fasta')
  in_data = read_ids( argsdict['if'] )
  out_dict = {}
  out_file = open(argsdict['oname']+"/step1_out.txt", 'w')
  
  all_seqs = {}
  for fasta in fasta_sequences:
    name, sequence = fasta.id, fasta.seq
    for element in in_data:
      if element in name[0:len(element)+3]:
        out_file.write( ">" + name + "\n" )
        out_file.write( str(sequence) + "\n" )
        all_seqs[name] = str(sequence)
  out_file.close()

  all_proteins = []
  req_proteins = read_ids( argsdict['pf'] )
  all_species = defaultdict(list )
  for key in all_seqs:
    pieces = key.split('_')
    all_proteins.append( pieces[-1] )
    all_species[ pieces[0] ].append( pieces[-1] )
  all_proteins = list( set( all_proteins))

  if argsdict['debug']:
    print "All proteins = ", all_proteins
    print "len(all proteins) = ", len(all_proteins)
    print "len(all seqs) = ", len(all_seqs)
    for key in all_species:
      print key, len(all_species[key])

  base_name = argsdict['oname']+'/tmp_'
  aln_files = []
  #for protein in all_proteins:
  for protein in req_proteins:
    f1s = base_name + protein + ".fa"
    f1a = base_name + protein + ".aln"
    aln_files.append( f1a )
    out_file = open(f1s , "w" )
    for key in all_seqs:
      pieces = key.split('_')
      if protein in pieces[-1]:
         out_file.write( ">" + key + "\n" )
         out_file.write( str(all_seqs[key]) + "\n" )
    if argsdict['debug']:
      print "Working on",key,pieces[-1]
      print "Working on",protein
    out_file.close()
    sysnulldev = open( os.devnull, 'w' )
    subprocess.call( ['muscle','-in',f1s,'-out',f1a], 
                     stdout=sysnulldev, stderr=sysnulldev)
    sysnulldev.close()

  aln_seqs = read_many_alns( aln_files )
  good_lens = nominal_lens(aln_seqs)
  f = open(argsdict['oname']+"/final_unfiltered.aln","w")
  for species in all_species:
    subseqs = query_species(aln_seqs, species)
    big_protein = ""
    for protein in req_proteins:
      goodseqs = query_protein(subseqs, protein)
      if len(goodseqs) > 1:
        print("Error: len(goodseqs)>1, that means the alignment failed somehow.")
        print("Error: Bailing out ...")
        exit(-1)
      elif len(goodseqs) < 1:
        print("Fixing missing sequence, adding dashes for species="+species+" and protein="+protein)
        big_protein += good_lens[ protein ] * '-'
      else:
        big_protein += str( goodseqs[ goodseqs.keys()[0] ] )
    f.write(">"+species+"\n")
    f.write(big_protein+"\n")
  f.close()
  if argsdict['debug']:
    print subseqs.keys()
    print good_lens

def filter_tree( argsdict ):
  full_aln_f = argsdict['oname'] + "/final_unfiltered.aln"
  filtered_aln_f = argsdict['oname'] + "/" + argsdict['oname'] + "_filtered_"
  filtered_aln_f+= str(round(argsdict['min_cov_position'],2))+".aln"

  seq_dict = read_many_alns( [ full_aln_f ] )
  full_aln_d = {}
  full_aa_d = {}
  for key in seq_dict:
    full_aln_d[key] = str( seq_dict[key] )
    #print key, len(str(seq_dict[key]))
  min_nt = int( argsdict['min_cov_position']*len(seq_dict) )
  print "at least these sequences per position to be considered ... ", min_nt
  clean_aln_d = {}
  mask_str = []
  for i in range( len(full_aln_d[key]) ):
    this_counter = 0
    for key in full_aln_d:
      this_sequence = full_aln_d[key]
      if this_sequence[i] != '-':
        this_counter += 1
    if this_counter >= min_nt:
      mask_str.append(1)
    else:
      mask_str.append(0)
  #
  f = open(argsdict['if']+"_mapping")
  map_d = {}
  for rawline in f:
    line = rawline.strip()
    pieces = line.split("\t")
    map_d[pieces[0]] = pieces[1]
  f.close()

  f = open(filtered_aln_f,'w')
  for key in full_aln_d:
    key2 = key.split("_")
    key2 = key2[0]
    #f.write(">"+key+"_"+map_d[key2]+"\n")
    f.write(">"+key+"\n")
    out_str = ""
    for i in range(len(mask_str)):
      if mask_str[i] == 1:
        out_str += full_aln_d[key][i]
    f.write(out_str+"\n")
  f.close()
  print len(mask_str), len(out_str)
  FastTreeCMD = 'FastTree -nt %s > %s_tree' % ( filtered_aln_f, filtered_aln_f )
  print("About to execute : "+FastTreeCMD)
  subprocess.call( FastTreeCMD, shell=True )
  #
  f = open( filtered_aln_f + "_tree" )
  xxx = Phylo.read( f, format='newick' )
  f.close()
  #Phylo.draw( xxx )

def main_seq( argsdict ):
  try:
    os.mkdir( argsdict['oname'] )
  except OSError:
    if argsdict['overwrite'] == '0':
      print "Error: Output directory already exists & overwrite is off, aborting."
      exit(-1)
    else:
      print "Notice: Output directory already exists & overwrite is on, continuing."

  get_species_from_db( argsdict )
  filter_tree( argsdict )

def main():

  parser = argparse.ArgumentParser(description='ML classifier try 1')

  parser.add_argument('--if', help='Input file (seq ids)', 
                      required=False, default='input.ids')
  parser.add_argument('--df', help='Ref database file (fasta)', 
                      required=False, default='phyeco.fa')
  parser.add_argument('--pf', help='Protein list file',
                      required=False, default='protein.ids')
  min_cov_position_help = """Minimum coverage for a given position of the gene 
                             alignment to be considered into the combined 
                             alignment, float val in [0,1]."""
  parser.add_argument('--min_cov_position', help = min_cov_position_help,
                      type = float, required=False, default=0.5)
  oname_help = """base name for output files and dir, it will create a dir 
                  called ONAME, and put intermediate files. Tree will be 
                  named ONAME.tree, final alignment will be called ONAME.aln"""
  parser.add_argument('--oname', 
                      help=oname_help, 
                      required=False, 
                      default="tree_out")
  overwrite_help = 'If output files exist, and overwrite is != 0, it will overwrite them'
  parser.add_argument('--overwrite', 
                      help=overwrite_help, 
                      required=False, default='0')
  parser.add_argument('--debug', help='Write debug info', 
                      required=False, type=bool, default=False)

  args = parser.parse_args()
  argsdict = vars(args)

  main_seq( argsdict )

if __name__ == "__main__":
  main()

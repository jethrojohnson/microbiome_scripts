'''
taxid2taxonomy.py - fetch full taxonomy for ncbi taxids

Purpose
_______

This scripts contains methods for retrieving full taxonomic information for one
or more NCBI taxids. 

Consider using Entez Direct for this: 
https://www.ncbi.nlm.nih.gov/books/NBK179288/

Using names.dmp and nodes.dmp
+++++++++++++++++++++++++++++

The NCBI taxonomy files names.dmp and nodes.dmp can be downloaded from
ftp.ncbi.nlm.nih.gov/pub/taxonomy/. If provided, then the requested taxonomic
levels will be retrieved from these files.

Please note that there are 44 unique taxonomic ranks specified in nodes.dmp, 
only generic ranks (see --levels defaults) can be returned using this script. 

If you get 'missing' kingdom, then consider extending the maximum number of 
recursive lookups allowed (you will have to edit the script for this). 

If you request taxonomic information for levels below the taxid provided, 
they will be returned as 'unclassified'


Accessing a BioSQL db
+++++++++++++++++++++

This script currently expects a MySQL version of BioSQL. 

NOT YET IMPLEMENTED


Usage:

cat taxids | taxid2taxonomy --names-dmp=/path/to/names.dmp --nodes.dmp=/path/to/nodes.dmp


'''


import os,sys,re
import collections
import cgatcore.experiment as E
import cgatcore.iotools as iotools


limit = 0
upper_limit=30

def fetchTaxonomy(out_dict, taxid, taxid2scname, taxid2parent, taxid2rank):
    '''Receive a taxid number and return k,p,c,o,g,s for that taxid'''
    for k,v in out_dict.items():
        out_dict[k] = 'unclassified'

    # HACK: To stop infinate recursion a finite limit is set
    global limit
    limit = 0
    
    def _fetchTax(taxid):
        global limit
        try:
            if taxid2rank[taxid].endswith('kingdom') or limit > upper_limit:
                if taxid2rank[taxid].endswith('kingdom'):
                    out_dict['kingdom'] = taxid2scname[taxid]
                return out_dict
            elif taxid2rank[taxid] in out_dict.keys():
                limit += 1
                out_dict[taxid2rank[taxid]] = taxid2scname[taxid]
                return _fetchTax(taxid2parent[taxid])
            else:
                limit += 1
                return _fetchTax(taxid2parent[taxid])
        except KeyError:
            # If value is missing, assume it must be first value
            for k in out_dict.keys():
                out_dict[k] = 'missing'
            return out_dict
    return _fetchTax(taxid)


def parseNCBINodes(names_dmp, nodes_dmp):
    '''
    Iterate over NCBI names and nodes files to create mapping dictionaries
    for recursive look up of taxids.
    '''
    # parse names.dmp (taxid | name | unique name | type of name)
    taxid2scname = {}
    for line in iotools.open_file(names_dmp):
        line = line.strip('\t|\n').split('\t|\t')
        if line[3] == 'scientific name':
            taxid2scname[line[0]] = line[1]
        else:
            continue

    E.info('taxid2scname is %fMb, with %i keys' \
            % (sys.getsizeof(taxid2scname)/1000000.0, len(taxid2scname.keys())))

    # parse nodes.dmp (taxid | parent taxid | rank |...)
    taxid2parent = {}
    taxid2rank = {}
    for line in iotools.open_file(nodes_dmp):
        line = line.strip('\t|\n').split('\t|\t')
        taxid2parent[line[0]] = line[1]
        taxid2rank[line[0]] = line[2]

    E.info('taxid2parent is %fMb, with %i keys' \
           % (sys.getsizeof(taxid2parent)/1000000.0, len(taxid2parent.keys())))
    E.info('taxid2rank is %fMb with %i keys' \
           % (sys.getsizeof(taxid2rank)/1000000.0, len(taxid2rank.keys())))

    return taxid2scname, taxid2parent, taxid2rank


def fetchBioSQLTaxonomy(ncbi_tax_id,
                        host='',
                        user='guest',
                        passwd='',
                        db='biosql'):
    '''NOT YET IMPLEMENTED 
    Fetch the full parent taxonomy (to Kingdom) for
       a given NCBI tax id number'''

    phylogeny = collections.OrderedDict()

    # connect to BioSQL database
    db = MySQLdb.connect(host=host, user=user, passwd=passwd, db=db)
    cur = db.cursor()

    # break if more than 50 iterations...
    safety = 50
    first = True
    node_rank = None
    while node_rank != 'superkingdom':
        safety -= 1
        if safety < 0:
            break

        # internal BioSQL taxon_ids don't correspond to ncbi taxon IDs
        # WARNING: There is overlap between the two. 
        field = 'taxon_id'
        if first:
            tax_id = ncbi_tax_id
            field = 'ncbi_taxon_id'
            first = False
            
        # fetch tax_id, parent_tax_id node_rank, scientific_name
        statement = ("SELECT "
                     "  taxon.taxon_id,"
                     "  taxon.parent_taxon_id,"
                     "  taxon.node_rank,"
                     "  taxon_name.name"
                     " FROM taxon JOIN taxon_name"
                     "  ON taxon.taxon_id = taxon_name.taxon_id"
                     " WHERE taxon.{} = {} "
                     "AND name_class = 'scientific name'".format(field, tax_id))
        cur.execute(statement)
        taxonomy = cur.fetchall()
        if len(taxonomy) == 1:
            old_tax_id, tax_id, node_rank, scientific_name = taxonomy[0]
            phylogeny[node_rank] = scientific_name
        elif len(taxonomy) == 0:
            E.warn('No taxonomy information for taxon ID %s beyond %s' \
                   % (ncbi_tax_id, node_rank))
            phylogeny['genus'] = 'unavailable'
            phylogeny['species'] = 'unavailable'
            
            return phylogeny

        else:
            raise ValueError('multiple entries for tax id %s: %s' \
                             % (ncbi_tax_id, str(taxonomy)))
            
    return phylogeny




def main(argv=None): 
    
    if not argv:
        argv = sys.argv

    parser = E.ArgumentParser(description=__doc__)
    
    parser.add_argument("--version", action="version", version="1.0")
    parser.add_argument("--names-dmp", dest="names_dmp", type=str,
                        help="Path to ncbi names.dmp file")
    parser.add_argument("--nodes-dmp", dest="nodes_dmp", type=str,
                        help="Path to ncbi nodes.dmp file")
    parser.add_argument("--header", action="store_true", dest="header",
                        help="Does input contain a header?")
    parser.add_argument("--missing", dest="missing", type=str,
                        help="Text to interpret as a missing taxid")
    parser.add_argument("--levels", dest="levels", type=str,
                        help="A comma-separated list of the taxonomic"\
                        " levels to be returned")
    parser.add_argument("--header-out", dest="header_out", action="store_true",
                        help="Include a header in output")

    parser.set_defaults(names_dmp=None,
                        nodes_dmp=None,
                        header=False,
                        header_out=False,
                        missing="missing",
                        levels="kingdom,phylum,class,order,family,genus,species,strain")

    (args) = E.start(parser, argv=argv)

    # Set dictionary for taxonomic levels to be returned
    t = collections.OrderedDict()
    for i in args.levels.split(','):
        t[i] = 'missing'

    # Output header if requested
    if args.header_out:
        args.stdout.write('\t'.join([x for x in args.levels.split(',')]) + '\n')

    # Recursive look up of names.dmp and nodes.dmp files
    if args.names_dmp or args.nodes_dmp: 
        assert args.names_dmp and args.nodes_dmp, \
            "Please supply path to both names.dmp and nodes.dmp files"

        E.info("Generating dictionaries for recursive taxonomic look up:")
        taxid2scname, taxid2parent, taxid2rank = parseNCBINodes(args.names_dmp,
                                                                args.nodes_dmp)

        header = args.header
        for tax_id in args.stdin:
            tax_id = tax_id.strip()
            if header:
                header = False
                continue
            
            if not tax_id == args.missing:
                t = fetchTaxonomy(t,
                                  tax_id,
                                  taxid2scname,
                                  taxid2parent,
                                  taxid2rank)

            args.stdout.write('\t'.join([t[x] for x in args.levels.split(',')])\
                              + '\n')

if __name__=="__main__":
    sys.exit(main(sys.argv))

import sys
import PSet
files = []
outfile = file( 'filesToProcess.txt', 'w')
for ifile in PSet.process.source.fileNames :    
    outfile.write('root://cms-xrd-global.cern.ch/' + ifile + '\n' )


outfile.close()

sys.argv.append('--files')
sys.argv.append('filesToProcess.txt')

print sys.argv

#from NtupleReader_fwlite import *

#NtupleReader_fwlite(sys.argv)

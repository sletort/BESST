'''
	Created on March 08, 2017

	@author: sletort

	This file is part of BESST.

	BESST is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	BESST is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with BESST.  If not, see <http://www.gnu.org/licenses/>.
'''

import unittest
import random


import sys
import os.path
from inspect import getsourcefile
cur_file = os.path.realpath( getsourcefile( lambda:0 ) )
sys.path.insert( 0, os.path.dirname( cur_file ) + '/../../BESST' )

# ---------------------------------------
def GenerateContig( name, pos, l_ctg_len_range ):
	direction = random.choice( [ True,False ] )
	length    = random.randrange( *l_ctg_len_range )
	return ( name, direction, pos, length, "" )

def GenerateScaffold( n, l_ctg_len_range, l_gap_stat ):
	l_scaff = []

	l_gaps  = [ int( random.gauss( *l_gap_stat ) ) for i in range( n-1 ) ]

	pos = 1
	for i in range( n-1 ):
		l_scaff.append( GenerateContig( i+1, pos, l_ctg_len_range ) )
		pos += l_scaff[-1][3] + l_gaps[i] # pos = pos + ctg_len + gap

	# dernier scaff
	l_scaff.append( GenerateContig( n, pos, l_ctg_len_range ) )

	return l_scaff

class Param( object ):
	def __init__( self, std_dev, max_overlap, info_file ):
		self.std_dev_ins_size   = std_dev
		self.max_contig_overlap = max_overlap
		self.information_file   = info_file

# ---------------------------------------

#~ class GenerateOutputTest( unittest.TestCase ):

	#~ def setUp( self ):
		#~ """Define the environnement mandatory for a test."""

	#~ def tearDown( self ):
		#~ """Destroy the environnement set in setUp method."""

	#~ def test_PrintOutput( self ):
		#~ """?"""

import GenerateOutput as GO

class ScaffoldTest_without_overlap( unittest.TestCase ):

	def setUp( self ):
		"""Build the list of contig tuples."""
		l_ctgs = [
			( 'a', True, 0, 10, 'a'*10 ),
			# tiny gap
			( 't', True, 12, 10, 't'*10 ),
			# empty gap
			( 't_inv', False, 22, 10, 't'*10 ),
			# big gap
			( 'c', True, 45, 10, 'c'*10 )
		]
		self.info_file = open( "infos", "w" )
		p = Param( 32, 10, self.info_file )
		self.scaff_name = "test1"
		self.o_scaff = GO.Scaffold( self.scaff_name, p, l_ctgs )

	def tearDown( self ):
		self.info_file.close()


	def test_isScaffold( self ):
		self.assertIsInstance( self.o_scaff, GO.Scaffold )

	def test_Fasta( self ):
		file = "test1.fasta"
		with open( file, 'w' ) as o_f:
			self.o_scaff.make_fasta_string( o_f )

		l_expected = [
			">" + self.scaff_name + "\n",
			'a'*10 + 'NN' + 't'*10 + 'n'
				+ 'a'*10 + 'N'*13 + 'c'*10 + "\n"
		]
		with open( file ) as f_in:
			l_done = f_in.readlines()
		self.assertListEqual( l_done, l_expected )

	def test_Agp( self ):
		file = "test1.agp"
		with open( file, 'w' ) as o_f:
			self.o_scaff.make_AGP_string( o_f )

		l_expected = [
			"\t".join([ str(x) for x in [ self.scaff_name,  1,10, 1, 'W', 'a', 1, 10, '+' ] ]) + "\n",
			"\t".join([ str(x) for x in [ self.scaff_name, 11,12, 2, 'N', 2, 'scaffold', 'yes', 'paired-ends' ] ]) + "\n",
			"\t".join([ str(x) for x in [ self.scaff_name, 13,22, 3, 'W', 't', 1, 10, '+' ] ]) + "\n",
			"\t".join([ str(x) for x in [ self.scaff_name, 23,23, 4, 'N', 1, 'scaffold', 'yes', 'paired-ends' ] ]) + "\n",
			"\t".join([ str(x) for x in [ self.scaff_name, 24,33, 5, 'W', 't_inv', 1, 10, '-' ] ]) + "\n",
			"\t".join([ str(x) for x in [ self.scaff_name, 34,46, 6, 'N', 13, 'scaffold', 'yes', 'paired-ends' ] ]) + "\n",
			"\t".join([ str(x) for x in [ self.scaff_name, 47,56, 7, 'W', 'c', 1, 10, '+' ] ]) + "\n"
		]
		with open( file ) as f_in:
			l_done = f_in.readlines()
		self.assertListEqual( l_done, l_expected )

	def test_Gff( self ):
		file = "test1.gff"
		with open( file, 'w' ) as o_f:
			self.o_scaff.make_GFF_string( o_f )

		l_expected = [
			"\t".join([ self.scaff_name, "besst_assembly", 'contig', '1','10', '.', '+', '.', "ID={0};Name={0}".format( 'a' ) ]) + "\n",
			"\t".join([ self.scaff_name, "besst_assembly", 'gap', '11','12', '.', '.', '.', "" ]) + "\n",
			"\t".join([ self.scaff_name, "besst_assembly", 'contig', '13','22', '.', '+', '.', "ID={0};Name={0}".format( 't' ) ]) + "\n",
			"\t".join([ self.scaff_name, "besst_assembly", 'gap', '23','23', '.', '.', '.', "" ]) + "\n",
			"\t".join([ self.scaff_name, "besst_assembly", 'contig', '24','33', '.', '-', '.', "ID={0};Name={0}".format( 't_inv' ) ]) + "\n",
			"\t".join([ self.scaff_name, "besst_assembly", 'gap', '34','46', '.', '.', '.', "" ]) + "\n",
			"\t".join([ self.scaff_name, "besst_assembly", 'contig', '47','56', '.', '+', '.', "ID={0};Name={0}".format( 'c' ) ]) + "\n"
		]
		with open( file ) as f_in:
			l_done = f_in.readlines()
		self.assertListEqual( l_done, l_expected )

class ScaffoldTest_with_overlap( unittest.TestCase ):

	def setUp( self ):
		"""Build the list of contig tuples."""
		l_ctgs = [
			( 'cplx1', False, 0, 40, 'a'*10+'c'*10+'g'*10+'t'*10 )
			# True overlap of 25 char
			#	gap is a little shorter
			,( 'cplx2', True, 16, 40, 'c'*5+'g'*10+'t'*10+'a'*10+'c'*5 ),
			# True overlap of 10 char
			#	gap is a little longer
		]
		self.info_file = open( "infos", "w" )
		p = Param( 32, 30, self.info_file )
		self.scaff_name = "test1"
		self.o_scaff = GO.Scaffold( self.scaff_name, p, l_ctgs )

	def tearDown( self ):
		self.info_file.close()

	def test_Fasta( self ):
		file = "test1.fasta"
		with open( file, 'w' ) as o_f:
			self.o_scaff.make_fasta_string( o_f )

		l_expected = [
			">" + self.scaff_name + "\n",
			'a'*10+'c'*10+'g'*10+'t'*10
				+ 'n'
				+ 'a'*10+'c'*5
				+ "\n"
		]
		with open( file ) as f_in:
			l_done = f_in.readlines()
		self.assertListEqual( l_done, l_expected )

	def test_Agp( self ):
		file = "test1.agp"
		with open( file, 'w' ) as o_f:
			self.o_scaff.make_AGP_string( o_f )

		l_expected = [
			  "\t".join([ str(x) for x in [ self.scaff_name,  1,40, 1, 'W', 'cplx1', 1, 40, '-' ] ]) + "\n"
			, "\t".join([ str(x) for x in [ self.scaff_name, 41,41, 2, 'N', 1, 'scaffold', 'yes', 'paired-ends' ] ]) + "\n"
			, "\t".join([ str(x) for x in [ self.scaff_name, 42,56, 3, 'W', 'cplx2', 26, 40, '+' ] ]) + "\n"
		]
		with open( file ) as f_in:
			l_done = f_in.readlines()
		self.assertListEqual( l_done, l_expected )

	def test_Gff( self ):
		file = "test1.gff"
		with open( file, 'w' ) as o_f:
			self.o_scaff.make_GFF_string( o_f )

		l_expected = [
			  "\t".join([ self.scaff_name, "besst_assembly", 'contig', '1','40', '.', '-', '.', "ID={0};Name={0}".format( 'cplx1' ) ]) + "\n"
			, "\t".join([ self.scaff_name, "besst_assembly", 'gap', '41','41', '.', '.', '.', "" ]) + "\n"
			, "\t".join([ self.scaff_name, "besst_assembly", 'contig', '42','56', '.', '+', '.', "ID={0};Name={0}".format( 'cplx2' ) ]) + "\n"
		]
		with open( file ) as f_in:
			l_done = f_in.readlines()
		self.assertListEqual( l_done, l_expected )
		#~ self.assertEqual( l_done[1], l_expected[1] )

if __name__ == '__main__':
	unittest.main()

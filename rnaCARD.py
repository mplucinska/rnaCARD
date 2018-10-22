
# wersja 10.01.2018
# funkcje: wyszukiwanie pasujacych spinek, zamykajacych motywy stemow oraz nakladajacch sie stemow

import math
import argparse
import subprocess
import os
import sys


class Arguments:
	def __init__(self):
		fasta = ""
		cond_1 = ""
		cond_2 = ""
		i = ""
		m = False
		n = False
		out = ""
		overlap_small = 0.6

	def parse(self):
		parser = argparse.ArgumentParser()
		parser.add_argument("-i", help = "input file", required = True, type = str)
		parser.add_argument("--mo", help = "minimum overlap", type = str, default = 0.6)
		#parser.add_argument("--cs", help = "closing stems", action='store_true', required = False)
		parser.add_argument("--os", help = "overlaping stems", action='store_true', required = False)
		parser.add_argument("--match", help = "find matching motifs", action='store_true', required = False)
		parser.add_argument("--mismatch", help = "find mismatching motifs", action='store_true', required = False)

		args = parser.parse_args()
		self.i = args.i
		self.overlap_small = args.mo
		#self.cs = args.cs
		self.os = args.os
		self.match = args.match
		self.mismatch = args.mismatch

class Structure:
	def _init_(self):
		self.shape = "" 	#structure in abstract shape notation (level 5)
		self.bracket = ""	#structure in dot-bracket notation
		self.shape_position = []	#start and end position for each element in shape [i] = [start, end]
		self.domains = []			#shape domains with loop sequences
		self.domains_position = []
		self.sequence = ""
		self.mismatched_domains = []	#domains which appears only in one structure 
		self.subdomain_shape_position = []
		self.subdomain_match = []
		self.mismatch = []
		self.match = []
		self.closed_motifs = []
		self.str_id = ''
	def find_bracket(self,p):
		c_o = 0
		c_c = 0
		x = "Z"
		if self.bracket[p] == "(":
			for i in range(p , len(self.bracket) + 1):
				if self.bracket[i] == "(":
					c_o += 1
				elif self.bracket[i] == ")":
					c_o = c_o - 1
				if c_o == 0:
					x = i
					break
		elif self.bracket[p] == ")":
			for i in list(reversed(range(0, p + 1))):
				if self.bracket[i] == ")":
					c_o += 1
				elif self.bracket[i] == "(":
					c_o -= 1
				if c_o == 0:
					x = i
					break
		return x
	def new_stem(self, a, b, a_p, b_p):
		lvl = 5
		sep = False
		for i in range(a, b + 1):
			if self.bracket[i] != ".":
				if self.bracket[i] != self.bracket[a]:
					sep = True
					break
		return sep
	def bracket_to_shape(self):
		self.shape = ""
		self.shape_position = {}
		st = 0
		beg = False 
		cl = False
		prev_i = 0
		lvl = 5

		for i in range(0, len(self.bracket)):
			if self.bracket[i] == "(":
				if beg == False:
					self.shape += "["
					st = i
					self.shape_position[len(self.shape) - 1] = [i]
					beg = True
					#print prev_pair
				elif beg == True:
					akt_pair = self.find_bracket(i)
					prev_pair = self.find_bracket(prev_i)
					if self.new_stem(akt_pair, prev_pair,prev_i, i) == True:
						beg = False
						self.shape_position[len(self.shape)-1].append(prev_i)
						self.shape += "["
						st = i
						self.shape_position[len(self.shape) - 1] = [i]
						beg = True
						prev_pair = self.find_bracket(i)
				if cl == True:
					cl = False
					self.shape_position[len(self.shape)-2].append(prev_i)

			elif self.bracket[i] == ")":
				if beg == True:
						#print i, prev_i
						beg = False
						self.shape_position[len(self.shape)-1].append(prev_i)
				if cl == False:
					self.shape += "]"
					self.shape_position[len(self.shape) - 1] = [i]
					cl = True
					akt_pair = self.find_bracket(i)
					prev_pair = self.find_bracket(prev_i)
				elif cl == True:
					akt_pair = self.find_bracket(i)
					prev_pair = self.find_bracket(prev_i)
					if self.new_stem(akt_pair, prev_pair, prev_i, i) == True:
						cl = False
						self.shape_position[len(self.shape)-1].append(prev_i)						
						self.shape += "]"
						self.shape_position[len(self.shape) - 1] = [i]
						cl = True
			#print self.bracket[i], i
			if self.bracket[i] !=  ".":
				prev_i = i
		try:
			self.shape_position[len(self.shape)-1].append(prev_i)
		except KeyError:
			return False
	def get_domains(self): #dividing shape into domains
		o = 0
		c = 0
		poz = 0
		o_p = 0
		c_p = 0
		loops = self.loop_seq()
		loop_num = 0
		self.domains_position =[]
		self.domains = []
		op = False
		for i, j in enumerate(self.shape):
			if j == "[" :
				if self.shape[i + 1] == "[":
					self.domains.append("[")
					self.domains_position.append(self.shape_position[i])
				else:
					op = True
			elif j == "]":
				if op == True:
					dom = "[" + loops[loop_num] + "]"
					self.domains.append(dom)
					self.domains_position.append([self.shape_position[i - 1][0], self.shape_position[i][1]])
					loop_num += 1
					op = False
				else:
					self.domains.append("]")
					self.domains_position.append(self.shape_position[i])
	def loop_seq(self):
		o = False
		op =0
		c = False
		cp = 0
		d = False
		p = 0
		sek_list = []
		for i in self.bracket:
			if i == '(':
				o = True
				op = p
			elif o == True and i == ')':
				c = True
				cp =p
			if o == True and c ==True:
				sek_list.append(self.sequence[op + 1 : cp])
				o = False
				c = False
			p += 1
		return sek_list

class Transcript:
	def _init_(self):
		self.id = ""
		self.sequence = ""
		self.domain_match = []
		self.mismatch_positions = []
		self.match_positions = []
		self.match_positions_small = []
		self.closed_match_position = []
		self.motifs_count = 0
		self.matched_motifs = {}
	def input_a(self):
		out_file = open("matched_motifs_out.txt", 'w')
		out_file_2 = open("matched_whole_transcripts.txt", 'w')
		out_file.close()
		out_file_2.close()

		with open(arg.i) as input_f:
			for i in input_f:
				if i.startswith(">"):
					self.id = i.strip()[1:]
					self.sequence = input_f.next().strip().split()[1]
					s1 = Structure()
					line = input_f.next().strip().split()	# structure 1
					s1.bracket = line[1]
					s1.str_id = line[0]
					s2 = Structure()
					line = input_f.next().strip().split()	# structure 2
					s2.bracket = line[1]
					s2.str_id = line[0]

					s1.sequence = self.sequence
					s2.sequence = self.sequence
					#print ">",self.id
					self.match_positions = []
					self.match_positions_small = []

					if s1.bracket_to_shape() != False and s2.bracket_to_shape() != False:
						if s1.bracket != s2.bracket:
							s1.get_domains() #returns s1.domains + s1.domains_position
							s2.get_domains()
							#matching motifs
							if arg.match:
								self.get_matched_hairpins_1(s1,s2) # pasujace hairpiny
								#self.get_matched_hairpins_2(s1,s2)
								#if arg.cs:
								#	self.get_matched_closing_stems(s1,s2) #stemy zamykajace dopasowane motywy
								if arg.os:
									self.get_matched_similar_stems(s1,s2) # podobne stemy
								self.print_match(s1,s2)
								self.matched_motifs_output(s1,s2)
								#self.print_mismatch(s1,s2)
							#mismatching motif
							elif arg.mismatch:
								self.get_matched_hairpins_1(s1,s2)
								#self.get_matched_closing_stems(s1,s2)
								if arg.os:
									self.get_matched_similar_stems(s1,s2)
								self.print_mismatch(s1,s2)
						#else:
						#	print "takie same!"
	
	def add_matching_motifs(self,s1,s2, pos_start, pos_end, shape_start, shape_end, shape2_start):
		#self.match_positions.append([pos_start, pos_end])
		#for k in range(pos_start, pos_end + 1):
		#	self.match_string[int(k)] = '1'

		#find small closed motifs
		dif = shape2_start - shape_start
		i = shape_start
		new_shape_start = i
		while i < shape_end + 1:
			if s1.domains[i] == "[":
				pair = self.find_pair(s1, i)
				if pair == shape_end:
					self.motifs_count += 1
					self.matched_motifs[self.motifs_count] = [[pos_start, pos_end], s1.domains[shape_start : shape_end + 1], [shape_start, shape_end], dif]
					self.list_of_opened_stems[i] = self.motifs_count
					i = shape_end + 1
				elif pair < shape_end:
					i = pair
				elif pair > shape_end:
					if  i == shape_start:
						self.motifs_count += 1
						self.matched_motifs[self.motifs_count] = [[pos_start, s1.domains_position[i][1]], s1.domains[new_shape_start : i + 1], [new_shape_start, i], dif]
						pos_start = s1.domains_position[i + 1][0]
						new_shape_start = i + 1
						self.list_of_opened_stems[i] = self.motifs_count
					elif i == new_shape_start:
						#print i
						self.motifs_count += 1
						self.matched_motifs[self.motifs_count] = [[pos_start, s1.domains_position[i][1]], s1.domains[new_shape_start : i + 1], [new_shape_start, i], dif]
						self.list_of_opened_stems[i] = self.motifs_count

						pos_start = s1.domains_position[i + 1][0]
						new_shape_start = i + 1
					elif i != new_shape_start:
						self.motifs_count += 1
						self.matched_motifs[self.motifs_count] = [[pos_start, s1.domains_position[i - 1][1]],s1.domains[new_shape_start : i], [new_shape_start, i], dif]
						pos_start = s1.domains_position[i][0]
						self.motifs_count += 1
						self.matched_motifs[self.motifs_count] = [[pos_start, s1.domains_position[i][1]], s1.domains[i : i + 1],[i,i], dif]
						pos_start = s1.domains_position[i + 1][0]
						new_shape_start = i + 1
						self.list_of_opened_stems[i] = self.motifs_count
					i += 1
			elif s1.domains[i] == "]":
				pair = self.find_pair(s1, i)
				if pair < shape_start:
					self.matched_motifs[self.list_of_opened_stems[pair]][0].append(s1.domains_position[i][0])
					self.matched_motifs[self.list_of_opened_stems[pair]][0].append(s1.domains_position[i][1])
					self.matched_motifs[self.list_of_opened_stems[pair]][1].append(" ]")
					self.matched_motifs[self.list_of_opened_stems[pair]][2].append(i)
					self.matched_motifs[self.list_of_opened_stems[pair]][2].append(i)
					d = [self.matched_motifs[self.list_of_opened_stems[pair]][3], dif]
					self.matched_motifs[self.list_of_opened_stems[pair]][3] = d
				else:
					self.motifs_count += 1
					self.matched_motifs[self.motifs_count] = [[pos_start, s1.domains_position[i][1]], s1.domains[new_shape_start : i + 1], [new_shape_start, i], dif]

				if i < shape_end:
					pos_start = s1.domains_position[i + 1][0]
					new_shape_start = i + 1
				i += 1
			else:
				i += 1

		if s1.domains[shape_end] != "[" and s1.domains[shape_end] != "]":
			self.motifs_count += 1
			self.matched_motifs[self.motifs_count] = [[pos_start, pos_end], s1.domains[new_shape_start : shape_end + 1], [new_shape_start, shape_end], dif]

	def matched_motifs_output(self, s1, s2):
		out_file = open("matched_motifs_out.txt", 'a')
		for i in self.matched_motifs:
			pos = self.matched_motifs[i][0]
			shape_pos = self.matched_motifs[i][2]
			shape_dif = self.matched_motifs[i][3]
			#print self.matched_motifs[i][1]
			shape = " ".join(self.matched_motifs[i][1])
			#print shape
			if len(pos) == 2:
				line_out = self.id + "\t" + str(i) + "\t" + str(pos[0]) + "\t" + str(pos[1]) + '\t' + shape + "\t" + self.sequence[pos[0] : pos[1] + 1] + "\t" + s1.bracket[pos[0] : pos[1] + 1] + "\t" + s2.bracket[pos[0] : pos[1] + 1] + "\t" + str(s1.domains_position[shape_pos[0]][0]) + "\t" + str(s1.domains_position[shape_pos[1]][1]) + "\t" + str(s2.domains_position[shape_pos[0] + shape_dif][0]) + "\t" + str(s2.domains_position[shape_pos[1] + shape_dif][1]) + "\t" + "\n"
				for k in range(pos[0], pos[1] + 1):
					self.match_string[int(k)] = str(i)
			else:
				line_out = self.id + "\t" + str(i) + "\t" + str(pos[0]) + " " + str(pos[2]) + "\t" + str(pos[1]) + " " + str(pos[3])  + '\t' + shape + "\t" + self.sequence[pos[0] : pos[1] + 1] + "&" + self.sequence[pos[2] : pos[3] + 1] + "\t" + s1.bracket[pos[0] : pos[1] + 1] + "&" + s1.bracket[pos[2] : pos[3] + 1] + "\t" + s2.bracket[pos[0] : pos[1] + 1] + "&" + s2.bracket[pos[2] : pos[3] + 1] + "\t" + str(s1.domains_position[shape_pos[0]][0]) + " " + str(s1.domains_position[shape_pos[2]][1]) + "\t" + str(s1.domains_position[shape_pos[1]][0]) + " " + str(s1.domains_position[shape_pos[3]][1]) + "\t" + str(s2.domains_position[shape_pos[0] + shape_dif[0]][0]) + " " + str(s2.domains_position[shape_pos[1] + shape_dif[0]][1]) + "\t" + str(s2.domains_position[shape_pos[2] + shape_dif[1]][0]) + " " +  str(s2.domains_position[shape_pos[3] + shape_dif[1]][1]) + "\n"
				for k in range(pos[0], pos[1] + 1):
					self.match_string[int(k)] = str(i)
				for k in range(pos[2], pos[3] + 1):
					self.match_string[int(k)] = str(i)
			out_file.write(line_out)
		out_file.close()
		
		out_file_2 = open("matched_whole_transcript.txt", 'a')
		out_file_2.write(self.id + "\t" + self.sequence + "\t" + s1.bracket + "\t" + s2.bracket + "\t" + " ".join(self.match_string) + "\n")

		'''	
		print ">", self.id
		print self.sequence
		print s1.bracket
		print s2.bracket
		print " ".join(self.match_string)
		'''

	def print_match(self, s1, s2):
		self.match_string = list("0" * len(self.sequence))
		start = False
		self.motifs_count = 0
		self.matched_motifs = {}
		self.list_of_opened_stems = {} #number od stem : number of motif
		for i,n in enumerate(s1.pairs):
			if n != 'x' and start == False:
				if s1.domains_position[i][0] == s2.domains_position[n][0]:
					start = True
					pos_start = s1.domains_position[i][0]
					shape_start = i
					shape2_start = n
					if s1.domains_position[i][1] == s2.domains_position[n][1]:
						pos_end = s1.domains_position[i][1]
						shape_end = i
					else:
						pos_end = s1.domains_position[i][1] if s1.domains_position[i][1] < s2.domains_position[n][1] else s2.domains_position[n][1]
						shape_end = i
						self.add_matching_motifs(s1, s2, pos_start, pos_end, shape_start, shape_end, shape2_start)
						start = False
				else:
					if s1.domains_position[i][1] == s2.domains_position[n][1]:
						start = True
						pos_start = s1.domains_position[i][0] if s1.domains_position[i][0] > s2.domains_position[n][0] else s2.domains_position[n][0]
						pos_end = s1.domains_position[i][1]
						shape_start = i
						shape2_start = n
						shape_end = i
					else:
						pos_start = s1.domains_position[i][0] if s1.domains_position[i][0] > s2.domains_position[n][0] else s2.domains_position[n][0]
						pos_end = s1.domains_position[i][1] if s1.domains_position[i][1] < s2.domains_position[n][1] else s2.domains_position[n][1]
						shape_start = i
						shape2_start = n
						shape_end = i
						self.add_matching_motifs(s1, s2, pos_start, pos_end, shape_start, shape_end, shape2_start)
						start = False
			elif start == True:
				if n == 'x':
					start = False
					self.add_matching_motifs(s1, s2, pos_start, pos_end, shape_start, shape_end, shape2_start)
				elif n != 'x':
					if n != s1.pairs[i - 1] + 1:
						start = True
						shape_end = i - 1
						pos_end = s1.domains_position[i - 1][1] if s1.domains_position[i - 1][1] < s2.domains_position[s1.pairs[i - 1]][1] else s2.domains_position[s1.pairs[i - 1]][1]
						self.add_matching_motifs(s1, s2, pos_start, pos_end, shape_start, shape_end, shape2_start)
						pos_start = s1.domains_position[i][0] if s1.domains_position[i][0] > s2.domains_position[n][0] else s2.domains_position[n][0]
						shape_start = i
						shape2_start = n
						shape_end = i
						pos_end = s1.domains_position[i][1] if s1.domains_position[i][1] < s2.domains_position[n][1] else s2.domains_position[n][1]
					#else:
					if s1.domains_position[i][0] == s2.domains_position[n][0]:
						start = True
						shape_end = i
						pos_end = s1.domains_position[i][1] if s1.domains_position[i][1] < s2.domains_position[n][1] else s2.domains_position[n][1]
					else:
						#dodac starty motyw
						#pos_end = s1.domains_position[i][1] if s1.domains_position[i][1] < s2.domains_position[n][1] else s2.domains_position[n][1]
						if i != shape_start:
							self.add_matching_motifs(s1, s2, pos_start, pos_end, shape_start, shape_end, shape2_start)
						shape_start = i
						shape2_start = n
						pos_start = s1.domains_position[i][0] if s1.domains_position[i][0] > s2.domains_position[n][0] else s2.domains_position[n][0]
						#dodac nowy motyw jesli koniec nie pasuje
						if s1.domains_position[i][1] == s2.domains_position[n][1]:
							start = True
							shape_end = i
							pos_end = s1.domains_position[i][1]
						else:
							shape_start = i
							shape2_start = n
							pos_start = s1.domains_position[i][0] if s1.domains_position[i][0] > s2.domains_position[n][0] else s2.domains_position[n][0]
							shape_end = i
							pos_end = s1.domains_position[i][1] if s1.domains_position[i][1] < s2.domains_position[n][1] else s2.domains_position[n][1]
							self.add_matching_motifs(s1, s2, pos_start, pos_end, shape_start, shape_end, shape2_start)
							start = False
		if start == True:
			self.add_matching_motifs(s1, s2, pos_start, pos_end, shape_start, shape_end, shape2_start)

	def print_mismatch(self, s1, s2):
		#print s1.pairs
		#print s2.pairs
		#print s1.domains
		#print s2.domains
		mismatch_string = list("0" * len(self.sequence))
		start = False
		self.mismatch_positions = []

		if s1.pairs[0] != 0 and s2.pairs[0] != 0:
			prev_match_1 = -1
			prev_match_2 = -1
			start = True
		for i,n in enumerate(s1.pairs):
			if n != 'x':
				if start == False:
					prev_match_1 = i
					prev_match_2 = n
					start = True
				elif start == True:
					if n > prev_match_2 + 1 and i == prev_match_1 + 1:
						#print prev_match_2 + 1, n - 1, "s2"
						pos_start = s2.domains_position[prev_match_2 + 1][0]
						pos_end = s2.domains_position[n - 1][1]
						self.mismatch_positions.append([pos_start, pos_end])
						for k in range(pos_start, pos_end + 1):
							mismatch_string[int(k)] = '1'			
					elif n == prev_match_2 + 1 and i > prev_match_1 + 1:
						#print prev_match_1 + 1, i - 1, "s1"
						pos_start = s1.domains_position[prev_match_1 + 1][0]
						pos_end = s1.domains_position[i - 1][1]
						self.mismatch_positions.append([pos_start, pos_end])
						for k in range(pos_start, pos_end + 1):
							mismatch_string[int(k)] = '1'
					elif n > prev_match_2 + 1 and i > prev_match_1 + 1:
						#print prev_match_1, prev_match_2, i - 1, n - 1
						pos_start = s1.domains_position[prev_match_1 + 1][0] if s1.domains_position[prev_match_1 + 1][0] < s2.domains_position[prev_match_2 + 1][0] else s2.domains_position[prev_match_2 + 1][0]
						pos_end = s1.domains_position[i - 1][1] if s1.domains_position[i - 1][1] > s2.domains_position[n - 1][1] else s2.domains_position[n - 1][1]
						self.mismatch_positions.append([pos_start, pos_end])
						for k in range(pos_start, pos_end + 1):
							mismatch_string[int(k)] = '1'
					prev_match_2 = n
					prev_match_1 = i
		if n == 'x':
			try:
				pos_start = s1.domains_position[prev_match_1 + 1][0] if s1.domains_position[prev_match_1 + 1][0] < s2.domains_position[prev_match_2 + 1][0] else s2.domains_position[prev_match_2 + 1][0]
			except IndexError:
				try:
					pos_start = s1.domains_position[prev_match_1][0] if s1.domains_position[prev_match_1][0] < s2.domains_position[prev_match_2 + 1][0] else s2.domains_position[prev_match_2 + 1][0]
				except IndexError:
					try:
						pos_start = s1.domains_position[prev_match_1 + 1][0] if s1.domains_position[prev_match_1 + 1][0] < s2.domains_position[prev_match_2][0] else s2.domains_position[prev_match_2][0]
					except IndexError:
						pos_start = s1.domains_position[prev_match_1][0] if s1.domains_position[prev_match_1][0] < s2.domains_position[prev_match_2][0] else s2.domains_position[prev_match_2][0]

			pos_end = s1.domains_position[-1][1] if s1.domains_position[- 1][1] > s2.domains_position[- 1][1] else s2.domains_position[- 1][1]
			self.mismatch_positions.append([pos_start, pos_end])
			for k in range(pos_start, pos_end + 1):
				mismatch_string[int(k)] = '1'
		'''
		print ">", self.id
		print self.sequence
		print s1.bracket
		print s2.bracket
		print " ".join(mismatch_string)	
		'''
	def stems_overlap(self, s1, s2, pos_1, pos_2): #chceck if pairs match or in range
		match = False
		for i in range(pos_1[0], pos_1[1]):
			pair_1 = s1.find_bracket(i)
			for j in range(pos_2[0], pos_2[1]):
				if i == j:
					pair_2 = s2.find_bracket(j)
					if pair_1 == pair_2 and pair_1 != 'Z':
						#if self.id[:-4] == 'ENST00000306589':
						#	print i, j, pair_1, pair_2
						match = True
						break
		return match
	def get_matched_closing_stems(self, s1, s2):
		added_new = True
		while added_new == True:
			start = False
			added_new = False
			for i, n in enumerate(s1.pairs):
				if n == 'x':
					if start == False and i < len(s1.pairs) - 1:
						if s1.pairs[i + 1] != "x":
							start_s1 = i
							start_s2 = s1.pairs[i + 1] - 1
							if s1.domains[start_s1] == "[" and s2.domains[start_s2] == "[":
								start = True
					elif start == True:
						if s1.pairs[i] == "x":
							#print i, s1.pairs[i - 1] + 1
							try:
								if s1.domains[i] == "]" and s2.domains[s1.pairs[i - 1] + 1] == "]":
									end_s1 = i
									end_s2 = s1.pairs[i - 1] + 1
									start = False
									if self.stems_overlap(s1,s2,start_s1, start_s2, end_1, end_2):
										s1.pairs[start_s1] = start_s2
										s2.pairs[start_s2] = start_s1
										s1.pairs[end_s1] = end_s2
										s2.pairs[end_s2] = end_s1
										added_new = True
									#print "".join(s1.domains[start_s1 : end_s1 + 1])
									#print "".join(s2.domains[start_s2 : end_s2 + 1])
								elif s1.domains[i] == "[" and s2.domains[s1.pairs[i - 1] + 1] == "[":
									if s1.pairs[i + 1] != "x":
										start_s1 = i
										start_s2 = s1.pairs[i - 1] + 1
										start = True
									else:
										start = False
								else:
									start = False
							except IndexError:
								pass
	def find_pair(self, s, n):
		op = 0
		cl = 0
		if s.domains[n] == "[":
			for i in range(n, len(s.domains)):
				if s.domains[i] == "[":
					op += 1
				elif s.domains[i] == "]":
					cl += 1
				if op == cl:
					return i
		else:
			for i in list(reversed(range(0, n + 1))):
				if s.domains[i] == "[":
					op += 1
				elif s.domains[i] == "]":
					cl += 1
				if op == cl:
					return i
	def get_matched_similar_stems(self, s1, s2):
		for i, n in enumerate(s1.pairs):
			if n == 'x':
				if s1.domains[i] == "[":
					pos_1 = s1.domains_position[i]
					for j, m in enumerate(s2.domains):
						pos_2 = s2.domains_position[j]
						#if self.id[:-4] == 'ENST00000312310':
						#	print pos_2, len(set(range(pos_1[0], pos_1[1])).intersection(range(pos_2[0], pos_2[1])))
						if m == "[" and len(set(range(pos_1[0], pos_1[1])).intersection(range(pos_2[0], pos_2[1]))) > 0:
							if self.stems_overlap(s1,s2,pos_1, pos_2):
								pair_1 = self.find_pair(s1, i)
								pair_2 = self.find_pair(s2, j)
								pos_1_cl = s1.domains_position[pair_1]
								pos_2_cl = s2.domains_position[pair_2]
								if len(set(range(pos_1_cl[0], pos_1_cl[1])).intersection(range(pos_2_cl[0], pos_2_cl[1]))) > 0:
									s1.pairs[i] = j
									s1.pairs[pair_1] = pair_2
									s2.pairs[j] = i
									s2.pairs[pair_2] = pair_1
	def get_matched_hairpins_1(self,s1,s2):
		s1.pairs = ['x']*len(s1.domains)
		s2.pairs = ['x']*len(s2.domains)
		#print s1.domains
		#print s2.domains
		for n, i in enumerate(s1.domains):
			if len(i) > 1:
				for m, j in enumerate(s2.domains):
					if len(j) > 1:
						if self.loops_overlap(s1,s2,i,j, s1.domains_position[n], s2.domains_position[m]):
							#print i, j
							#print s2.bracket[s2.domains_position[m][0] : s2.domains_position[m][1] + 1]
							s1.domains[n] = s2.domains[m]
							s1.pairs[n] = m
							s2.pairs[m] = n
		#print s1.pairs
		#print s2.pairs
	def get_matched_hairpins_2(self,s1,s2):
		s1.pairs = ['x']*len(s1.domains)
		s2.pairs = ['x']*len(s2.domains)
		#print s1.domains
		#print s2.domains
		for n, i in enumerate(s1.domains):
			if len(i) > 1:
				for m, j in enumerate(s2.domains):
					if i == j :
						#print s2.bracket[s2.domains_position[m][0] : s2.domains_position[m][1] + 1]
						s1.pairs[n] = m
						s2.pairs[m] = n
		#print s1.pairs
		#print s2.pairs
	def loops_overlap(self, s1, s2, a, b, a_pos, b_pos):
		if a == b and a_pos == b_pos:
			#if len(set(range(a_pos[0], a_pos[1])).intersection(range(b_pos[0], b_pos[1]))) > 0:
			return True
		else:
			start_1 = "A"
			end_1 = "A"
			start_2 = "A"
			end_2 = "A"
			for j in range(a_pos[0], a_pos[1] + 1):
				if s1.bracket[j] == "(":
					start_1 = j + 1 
				elif s1.bracket[j] == ")":
					end_1 = j - 1
					break
			for j in range(b_pos[0], b_pos[1] + 1):
				if s2.bracket[j] == "(":
					start_2 = j + 1 
				elif s2.bracket[j] == ")":
					end_2 = j - 1
					break
			#if self.id[:-4] == 'ENST00000306589':
			#	print a, b, a_pos, b_pos, start_1, end_1, start_2, end_2
			if start_1 != "A" and start_2 != "A" and  end_1 != "A" and end_2 != "A":
				if end_1 - start_1 < end_2 - start_2:
					small_loop = end_1 - start_1
				else:
					small_loop = end_2 - start_2
				min_overlap = small_loop * arg.overlap_small
				if start_2 <= end_1 and start_2 >= start_1 and end_2 >= end_1:
					if end_1 - start_2 + 1 >= min_overlap + 1:
						return True
				elif start_1 <= end_2 and start_1 >= start_2 and end_1 >= end_2:
					if end_2 - start_1 + 1 >= min_overlap + 1:
						return True
				elif start_2 >= start_1 and end_2 <= end_1:
					return True
				elif start_1 >= start_2 and end_1 <= end_2:
					return True
			else:
				return False

arg = Arguments()
arg.parse()
Transcript().input_a()
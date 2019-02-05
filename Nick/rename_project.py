#------------------------------------------------------------------------------
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#------------------------------------------------------------------------------
# Filename: FindAndReplace.py
# Purpose:  Simple find and replace string in files (recursive) script
# Usage:	python FindAndReplace.py [Old String] [New String]
#				[File Filters(ex/default:".txt,.html,.erb")] [Directory To Check]
# Requirement: Files must be text (non-binary) files
#              (this is why we force you to pick the file pattern/filter)
# WARNING: This will overwrite files matching [File Filters]. All occurrences of [Old String]
#			will be replaced with [New String]. Make sure you really, really want to do this.
#------------------------------------------------------------------------------

import os
import sys
import traceback

def usage():
	print('Usage: python FindAndReplace.py [Old String] [New String] ' \
		  '[File Filters(default:".txt,.html,.erb")] [Directory To Check(.)]')

def replaceStringInFile(fileName, oldStringToReplace, newString):

	if not(os.path.isfile(fileName) and os.access(fileName, os.W_OK)):
		print("WARNING: Skipping...File does not exist or and is not writeable:" + fileName)
		return False

	fileUpdated = False

	# credit/taken/adapted from: http://stackoverflow.com/a/4128194
	# Read in old file
	with open(fileName, 'r', errors='ignore') as f:
		newlines = []
		for line in f.readlines():
			if (oldStringToReplace in line) :
				fileUpdated = True
			line = line.replace(oldStringToReplace, newString)
			newlines.append(line)

	# Write changes to same file
	if fileUpdated :
		print("String Found and Updating File: " + fileName)
		try:
			with open(fileName, 'w') as f:
				for line in newlines:
					f.write(line)
		except:
			print('ERROR: Cannot open/access existing file for writing: ' + fileName)

	return fileUpdated


def rename_filenames(workingDIR, old_name, new_name):
	counter = 0
	path = workingDIR
	for file in os.listdir(path):
		print(file, old_name, new_name)
		if old_name in file:
	            counter += 1
	            os.rename(path + "/"+file, path + "/"+file.replace(old_name, new_name))
	if counter == 0:
	    print("No file has been found")


def main():

	try:

		if len(sys.argv) < 4:
			usage()
			# old/new string required parameters, exit if not supplied
			sys.exit(-1)
		else:
			oldString = sys.argv[1]
			newString = sys.argv[2]
			path = sys.argv[3]+"/"+oldString

		print('[Old String]         : ' + oldString)
		print('[New String]         : ' + newString)
		print('[Directory To Check] : ' + path)

		if not os.path.exists(path):
			raise Exception("Selected path does not exist: " + path)

		# Walks through directory structure looking for files matching patterns
		#matchingFileList = \
		#	[os.path.join(dp, f) \
		#		for dp, dn, filenames in os.walk(path) #\
		#			for f in filenames \
		#				if os.path.splitext(f)[1] in patterns]

		matchingFileList = list()
		skip_folders=['c-sstar/','16s/', 'FASTQs/', 'kraken/', 'plasmid_on_plasmidAssembly/', 'removedAdapters/', 'rst2/', '/trimmed/', '.fastq']

		for (dirpath, dirnames, filenames) in os.walk(path):
			matchingFileList += [os.path.join(dirpath, file) for file in filenames]

		filesSkipped = 0
		for counter in range(0, len(matchingFileList)):
			for skip in skip_folders:
				if skip in matchingFileList[counter]:
					matchingFileList[counter] = ""
					filesSkipped+=1
					break

		matchingFileList = list(filter(None, matchingFileList))

		print('Files found: ' + str(len(matchingFileList)))
		print('Files skipped: ' + str(filesSkipped))
		fileCount = 0
		filesReplaced = 0

		for currentFile in matchingFileList:
			fileCount+=1
		#	print("Looking at", currentFile)
		#	for skip in skip_folders:
		#		print("Testing skip:", skip)
		#		if skip in strcurrentFile:
		#			print("Skipped", strcurrentFile)
		#			filesSkipped+=1
		#			break
		#		else:
			fileReplaced = replaceStringInFile(currentFile, oldString, newString)
			if fileReplaced:
				filesReplaced+=1

		print("Total Files Searched         : " + str(fileCount))
		print("Total Files Skipped          : " + str(filesSkipped))
		print("Total Files Replaced/Updated : " + str(filesReplaced))

		rename_filenames(path, oldString, newString)

		os.rename(path, sys.argv[3]+"/"+newString)

	except Exception as err:
		print(traceback.format_exception_only(type(err), err)[0].rstrip())
		sys.exit(-1)

if __name__ == '__main__':
	main()

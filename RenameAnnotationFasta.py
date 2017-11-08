import os
import string
import sys

"""
>NODE_1_length_1203982_cov_41.0681_1 # 578 # 754 # 1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=AGGAG;rbs_spacer=5-10bp;gc_cont=0.266
MKKIIIGVFFGILILVIVGLYVSYVYMNQSFIWVAGILLLLGTLFNWFLYNKYLSKKT*
>NODE_1_length_1203982_cov_41.0681_2 # 805 # 1422 # 1 # ID=1_2;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.358
MLMSLSVLPLATYASETENTPTESYGGENFIATQTGNTLVIEDKKTGETVKIEMNDEENG
VITSDDGTIENVHRDEEGNVYVDNELELEAPPLDIEDGINIATQPRLLKASKWIYVQTTK
YNTTTQGNMRSLALGILSFMPITGPIFGIVAIIDAARSMGAKTLYVRVKQYRTSGYQFYK
YDSYYYANASLTKLVKKTSQTKRMW*
>NODE_1_length_1203982_cov_41.0681_3 # 1732 # 2076 # 1 # ID=1_3;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.287
MQKDTYSGIRLSVIQLIDSKKENLKENNIKLTIIKDEKDGYVVELDNDKCMAEIVVEEPT
YAPYRYISFEVVSLMDGKVKIIYSWYDDETSQWSDIEKELNKGIQFLNNFKTME*
>NODE_1_length_1203982_cov_41.0681_4 # 2483 # 2725 # 1 # ID=1_4;partial=00;start_type=GTG;rbs_motif=AGGAG;rbs_spacer=5-10bp;gc_cont=0.292
MKRIILTSTLIVWTIVCIYMSISMVSNNTGIAFPIWLHIILLICFLATGIVNVKKKEYLW
STMLFEGVLVVLLSLIIVLV*
"""

infilename = sys.argv[1]
tmpfilename = infilename+'.tmp'
logfilename = infilename+'.stat'
samplename = infilename.split('_spades_')[0] # attenzione qui

proteincounter = 1

with open(infilename) as infile, open(tmpfilename, 'wb') as tmpfile, open(logfilename, 'wb') as statfile:
	for fasta in infile.read().split(">")[1:]:
		title = fasta.split("\n", 1)[0]
		sequence = fasta.split("\n", 1)[1]
		tmpfile.write(">"+samplename+'_prot'+str(proteincounter)+'\n'+sequence.strip()+'\n')
		statfile.write(samplename+'_prot'+str(proteincounter)+'\t'+title.strip()+'\n')
		proteincounter += 1
        
os.system("mv %s %s" % (tmpfilename, infilename))

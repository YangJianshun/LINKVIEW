#!/usr/bin/env python
## _*_coding:utf-8_*_
from __future__ import division
import argparse
import re
from fractions import Fraction
import sys
import os
binpath=os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(binpath)
import interval

parser = argparse.ArgumentParser()
parser.add_argument("input",help='The input file which contains location relationships. Can be a text file, format like: "chr1 start1 end1 chr2 start2 end2 color:opacity". It can also be a blast output file.')
#parser.add_argument('type',help='Type of input file, 0 for link.txt format like "chr1 start1 end1 chr2 start2 end2",1 for output file of blast format as "query subject %identity alignment_length mismatches gap_opens q.start q.end s.start s.end evalue bit score"')
parser.add_argument('-t','--type',help='Type of input file, 0 for link.txt format like "chr1 start1 end1 chr2 start2 end2 color:opacity", 1 for output file of blast format as "query subject identity alignment_length mismatches gap_opens q.start q.end s.start s.end evalue bit score"\tdefault=0')
parser.add_argument('-hl','--highlight',help='highlight.txt, Highlight section, format like chr start end [color]')
parser.add_argument('-k','--karyotype',help="Karyotype.txt, If you don't specify this file, it will automatically sort and display.format like chr1:[start1:end1] chr2...")
parser.add_argument('--opacity',help='opacity, default=0.5')
parser.add_argument('--color',help='color, default=green')
parser.add_argument('--svg_height',help='height of svg, default=800')
parser.add_argument('--svg_width',help='width of svg, default=1200')
parser.add_argument('--svg_thickness',help='thickness of chromosome, default=15')
parser.add_argument('-n','--svg_no_label',action="store_true",help='Do not show labels')
parser.add_argument('--svg_font_size',help='font size of the label, default=18')
parser.add_argument('--svg_label_angle',help='label rotation angle, default=0')
parser.add_argument('--svg_space',help='The proportion of white space left and right, default=0.2')
parser.add_argument('-s','--show_pos_with_label',action="store_true",help='Display location information on the label')
parser.add_argument('--scale',help='example:5k default=auto')
parser.add_argument('-o','--output',help='output file prefix, default=linkview_output')
parser.add_argument('--min_identity',help='if your input file is a output file of blast, the min_identy to filter the input file, default=95')
parser.add_argument('--min_alignment_length',help='if your input file is a output file of blast, the min_alignment_length to filter the input file, default=200')
parser.add_argument('--max_evalue',help='if your input file is a output file of blast, the max_evalue to filter the input file, default=1e-5')
parser.add_argument('--min_bit_score',help='if your input file is a output file of blast, the min_bit_score to filter the input file, default=5000')
parser.add_argument('--Gap_length',help='Length of gap between two segments per line,if > 1,It represents Physical length, if<1,It represents total_length_of_this_line * this. default=0.2')
parser.add_argument('--chr_len',help='chr_len.txt, format as: chr1 150000')
parser.add_argument('-P','--Parameter',help='Specify the parameters for each row separately in a file, if needed. E.g svg_font_size=20 show_pos_with_label=0 svg_no_label=1 Gap_length=500')


args=parser.parse_args()
if not args.type: args.type=0
else:args.type=int(args.type)
if not args.opacity: args.opacity=0.5
if not args.color: args.color='green'
if not args.svg_height: args.svg_height=800
else: args.svg_height = int(args.svg_height)
if not args.svg_width: args.svg_width=1200
else: args.svg_width= int(args.svg_width)
if not args.svg_thickness: args.svg_thickness=15
else: args.svg_thickness=int(args.svg_thickness)
if not args.svg_font_size: args.svg_font_size=18
if args.svg_label_angle:args.svg_label_angle=360-float(args.svg_label_angle)

if not args.svg_space: args.svg_space=0.2
else: args.svg_space=float(args.svg_space)
if not args.output: args.output='linkview_output'

if not args.Gap_length: args.Gap_length=0.2
if not args.min_identity: args.min_identity=95
else: args.min_identity=float(args.min_identity)
if not args.min_alignment_length: args.min_alignment_length=200
else: args.min_alignment_length=int(args.min_alignment_length)
if not args.max_evalue: args.max_evalue=1e-5
else: args.max_evalue=float(args.max_evalue)
if not args.min_bit_score: args.min_bit_score=5000
else: args.min_bit_score=float(args.min_bit_score)
relations=[]
posofchro={}
f=open(args.input,'r')
for line in f:
    if line.startswith('#'):continue
    line=line.strip()
    opacity=''
    color=''
    if(args.type==0):
        chr1,start1,end1,chr2,start2,end2,color = re.match(r'(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*(\S+)?',line).group(1,2,3,4,5,6,7)
        if color and ':' in color: color,opacity=color.split(':')
    elif(args.type==1):
        arr=line.split('\t')
        identity,alignment_length,evalue,bit_score=float(arr[2]),int(arr[3]),float(arr[10]),float(arr[11])
        if(not (identity>args.min_identity and alignment_length>args.min_alignment_length and evalue<args.max_evalue and bit_score>args.min_bit_score ) ):continue
        chr1,start1,end1,chr2,start2,end2=arr[0],arr[6],arr[7],arr[1],arr[8],arr[9]
        
    else:print('[Error] Unkonwn input type {}, only be 0 or 1'.format(args.type))
    if not opacity:opacity=args.opacity
    if not color:color=args.color
    start1=int(start1);end1=int(end1);start2=int(start2);end2=int(end2);
    posofchro.setdefault(chr1,[]);posofchro.setdefault(chr2,[]);
    posofchro[chr1]+=[start1,end1];posofchro[chr2]+=[start2,end2];
    relations.append({'chr1':chr1,'start1':start1,'end1':end1,'chr2':chr2,'start2':start2,'end2':end2,'color':color,'opacity':opacity})
if(args.chr_len and os.path.exists(args.chr_len)):
    chr_len={}
    f=open(args.chr_len,'r')
    for line in f:
        line=line.strip()
        chr_name,length=re.match(r'^(\S+)\s+(\S+)\s*$',line).group(1,2)
        chr_len[chr_name]=int(length)
    f.close()


class chro():
    def __init__(self,name,start,end,is_full):
        self.name=name
        self.start=start
        self.end=end
        self.length=end-start+1
        self.is_full=is_full
        self.left=None
        self.right=None
        self.top=None
        self.level=None
    def coordinate(self,x,is_up,is_start=False):
        if self.left and self.right and self.top:
            if is_up: y=self.top
            else: y=self.top+args.svg_thickness
            if is_start:x=self.left+(x-self.start)*scale
            else: x=self.left+(x-self.start+1)*scale
            return (x,y)
        else: return None

chros={}
for c in posofchro:
    chros[c]=chro(c,min(posofchro[c]),max(posofchro[c]),[True,True])
order={}
auto_start_end=[];
if(args.karyotype and os.path.exists(args.karyotype)):
    level = 0
    k=open(args.karyotype,'r')
    for line in k:
        if line.startswith('#'): continue
        order.setdefault(level,[])
        line=line.strip()
        chrs=re.split(r'\s+',line)
        for c in chrs:
            arr=c.split(':')
            name=arr[0]
            if(len(arr)==3 and re.match(r'\d+',arr[1]) and re.match(r'\d+',arr[2]) and int(arr[1])<int(arr[2])):
                chros[name]=chro(name,int(arr[1]),int(arr[2]),[True,True])
            else:auto_start_end.append(name)
            order[level].append(name)
        level+=1

    for i in range(len(relations)):
        relation=relations[i]
        relations[i].setdefault('display',True)
        chr1=relation['chr1'];chr2=relation['chr2']
        #print(relation)
        if not any([chr1 in x for x in order.values()]):relations[i]['display'] = False;continue
        if not any([chr2 in x for x in order.values()]):relations[i]['display'] = False;continue
        #if interval.relation([chros[chr1].start,chros[chr1].end],[relation['start1'],relation['end1']])==1: relations[i]['display'] = False;continue
        #if interval.relation([chros[chr2].start,chros[chr2].end],[relation['start2'],relation['end2']])==1: relations[i]['display'] = False;continue
        reverse = False
        if(relation['start1']-relation['end1'])*(relation['start2']-relation['end2'])<0: reverse=True;
        tmp_its1=interval.intersection([relation['start1'],relation['end1']] , [chros[relation['chr1']].start,chros[relation['chr1']].end])
        if not tmp_its1: relations[i]['display'] = False;continue;
        tmp_start1,tmp_end1 = sorted([relation['end1'],relation['start1']])
        len1 = tmp_end1-tmp_start1+1
        tmp_its1=[Fraction((tmp_its1[0]-tmp_start1),len1),Fraction((tmp_its1[1]-tmp_start1+1),len1)]
        tmp_its2=interval.intersection([relation['start2'],relation['end2']] , [chros[relation['chr2']].start,chros[relation['chr2']].end])
        #print('tmp_its2',tmp_its2)
        if not tmp_its2: relations[i]['display']=False;continue;
        tmp_start2,tmp_end2 = sorted([relation['end2'],relation['start2']])
        len2 = tmp_end2-tmp_start2+1
        tmp_its2=[Fraction((tmp_its2[0]-tmp_start2),len2),Fraction((tmp_its2[1]-tmp_start2+1),len2)]
        #print('tmp_its2',tmp_its2)
        #print(tmp_its1,tmp_its2)
        if(reverse == False):
            tmp_its1=tmp_its2=interval.intersection(tmp_its1,tmp_its2)
            if not tmp_its1: relations[i]['display'] = False;continue;
        else:
            tmp_its1=[1-tmp_its1[1],1-tmp_its1[0]]
            tmp_its2=interval.intersection(tmp_its1,tmp_its2)
            if not tmp_its2: relations[i]['display'] = False;continue;
            tmp_its1=[1-tmp_its2[1],1-tmp_its2[0]]
            tmp_start1=min(tmp_start1,tmp_end1)
            tmp_start2=min(tmp_start2,tmp_end2)
            tmp_its1.reverse()
        relations[i]['start1'] = int(tmp_start1 + tmp_its1[0]*len1)
        relations[i]['start2'] = int(tmp_start2 + tmp_its2[0]*len2)
        relations[i]['end1'] = int(tmp_start1 + tmp_its1[1]*len1 -1)
        relations[i]['end2'] = int(tmp_start2 + tmp_its2[1]*len2 -1)
    relations=[x for x in relations if x['display']==True]
    for c in auto_start_end:
        posofchro[c]=[]
    for relation in relations:
        if relation['chr1'] in auto_start_end:
            posofchro[relation['chr1']].append(relation['start1'])
            posofchro[relation['chr1']].append(relation['end1'])
        if relation['chr2'] in auto_start_end:
            posofchro[relation['chr2']].append(relation['start2'])
            posofchro[relation['chr2']].append(relation['end2'])
    for c in auto_start_end:
        chros[c].start=min(posofchro[c])
        chros[c]=chro(c,min(posofchro[c]),max(posofchro[c]),[True,True])
    #print('-------')
else:
#====Arrange the order of chros on the image. If you specify a Karyotype file, this step is omitted.===
    link={}
    for relation in relations:
        chr1=relation['chr1'];chr2=relation['chr2']
        link.setdefault(chr1,[]);link.setdefault(chr2,[])
        if not chr2 in link[chr1]:link[chr1].append(chr2)
        if not chr1 in link[chr2]:link[chr2].append(chr1)


    chro_ordered=[]
    first=min(link.keys(),key=lambda x:len(link[x]))
    chro_ordered.append(first)
    order[0]=[first]
    def find_next(pos,*chrs):
        global chro_ordered
        nexts=[]
        for chro in chrs:
            for next in link[chro]:
                if not next in chro_ordered:
                    nexts.append(next)
            chro_ordered += nexts
            if nexts: order[pos]=nexts
            find_next(pos+1,*nexts)

    find_next(1,first)

#======================================================================================================
#每一行分别的参数
svg_font_size={}
show_pos_with_label={}
svg_no_label={}
svg_label_angle={}
Gap_length={}

max_len=0
lens={}
middle_spaces={}
for line,chrs in order.items():
    num=len(chrs)
    length=0
    for chr in chrs:length+=chros[chr].length;
    svg_font_size[line]=args.svg_font_size
    show_pos_with_label[line]=args.show_pos_with_label
    svg_no_label[line]=args.svg_no_label
    svg_label_angle[line]=args.svg_label_angle
    Gap_length[line]=args.Gap_length
    lens[line]=length
if(args.Parameter and os.path.exists(args.Parameter)):
    P=open(args.Parameter,'r')
    level=0
    for line in P:
        m1=re.search(r'svg_font_size=(\d+(\.\d+)?)',line)
        if m1:svg_font_size[level]=m1.group(1)
        m2=re.search(r'show_pos_with_label=(\d+(\.\d+)?)',line)
        if m2:show_pos_with_label[level]=m2.group(1)
        m3=re.search(r'svg_no_label=(\d+(\.\d+)?)',line)
        if m3:svg_no_label[level]=m3.group(1)
        m4=re.search(r'Gap_length=(\d+(\.\d+)?)',line)
        if m4:
            Gap_length[level]=float(m4.group(1))
        
        m5=re.search(r'svg_label_angle=(\d+(\.\d+)?)',line)
        if m5:svg_label_angle[level]=360-float(m5.group(1))
        level+=1
for level,chrs in order.items(): #对middle_spaces[level]进行赋值
    if(Gap_length[level]>1):gap_len=int(Gap_length[level])
    else:gap_len=lens[level]*float(Gap_length[level])
    middle_spaces[level]=gap_len
    lens[level]=lens[level]+gap_len*(len(order[level])-1)

max_len=max(lens.values())
if args.svg_label_angle and len(svg_label_angle)>1:svg_label_angle[max(svg_label_angle.keys())]=360-svg_label_angle[max(svg_label_angle.keys())]

space_single=(args.svg_width)*args.svg_space/2
scale=args.svg_width*(1-args.svg_space)/max_len
L=(max_len/5)
if(args.scale):
    L,suf=re.match(r'(\d+)([km]?)',args.scale).group(1,2)
    L=int(L)
    if suf=='k':L=L*1000
    elif suf=='m':L=L*1000000

def get_L_str(L):
    i=10
    while 1:
         result=int(L/i)
         if not result: break
         n,unit=result,i
         i*=10
    L1=n*unit
    if(unit>=10000000):suffix='0 Mb'
    elif(unit>=1000000):suffix=' Mb'
    elif(unit>=100000):suffix='00 Kb'
    elif(unit>=10000):suffix='0 Kb'
    elif(unit>=1000):suffix=' Kb'
    elif(unit>=100):suffix='00 bp'
    elif(unit>=10):suffix='0 bp'
    else:suffix=' bp'
    s=str(n)+suffix
    return(L1,s)

L,S=get_L_str(L)
space_vertical=args.svg_height/(len(order)+1)
top=0
line_length=args.svg_width*(1-args.svg_space)/40
for line,chrs in order.items():
    left_start=(args.svg_width-lens[line]*scale)/2
    top+=space_vertical
    for c in chrs:
        if(args.chr_len and os.path.exists(args.chr_len) and chros[c].name in chr_len):
            chros[c].is_full=[chros[c].start==1,chros[c].end==chr_len[chros[c].name]]
        chros[c].level=line
        chros[c].top=top
        chros[c].left=left_start
        chros[c].right=left_start+chros[c].length*scale
        if not chros[c].is_full[0]: chros[c].left+=line_length;chros[c].right+=line_length
        left_start=chros[c].right+middle_spaces[line]*scale
        if not chros[c].is_full[1]: left_start+=line_length

#==================creat svg=================
max_level=len(order)-1
svg='<svg width="{}" height="{}" xmlns="http://www.w3.org/2000/svg" version="1.1">\n'.format(args.svg_width,args.svg_height)
for c in chros:
    if not any([c in x for x in order.values()]):continue
    angle_text=''
    if not chros[c].is_full[0]:
        svg+='''<line x1="{}" y1="{}" x2="{}" y2="{}" style="stroke:black;stroke-width:3"/>'''.format(chros[c].left-line_length,chros[c].top+args.svg_thickness/2,chros[c].left,chros[c].top+args.svg_thickness/2)
    svg+='''<rect x="{}" y="{}" width="{}" height="{}" fill="grey" />\n'''.format(chros[c].left,chros[c].top,chros[c].right-chros[c].left,args.svg_thickness)
    if not chros[c].is_full[1]:
        svg+='''<line x1="{}" y1="{}" x2="{}" y2="{}" style="stroke:black;stroke-width:3"/>'''.format(chros[c].right,chros[c].top+args.svg_thickness/2,chros[c].right+line_length,chros[c].top+args.svg_thickness/2)

    if(chros[c].level==0):label_x=chros[c].left;label_y=chros[c].top-args.svg_thickness/2
    elif(chros[c].level==max_level):
        label_x=chros[c].left;label_y=chros[c].top+args.svg_thickness*2;
        if(svg_label_angle[chros[c].level]):pass
    else:label_x=chros[c].right+8;label_y=chros[c].top+args.svg_thickness
    if not int(svg_no_label[chros[c].level]):
        label=chros[c].name
        if int(show_pos_with_label[chros[c].level]):label=label+'\t'+'{:,}'.format(chros[c].start)+' ~ '+ '{:,}'.format(chros[c].end)+' bp'
        angle_text=' transform="rotate({},{} {})"'.format(svg_label_angle[chros[c].level],label_x,label_y)
        svg+='''<text x="{}" y="{}" fill="black"  font-size="{}"{}>{}</text>\n'''.format(label_x,label_y,svg_font_size[chros[c].level],angle_text,label)

for relation in relations:
    chr1=relation['chr1'];chr2=relation['chr2']
    is_up=chros[chr1].level>chros[chr2].level
    vertex1 = chros[chr1].coordinate(relation['start1'],is_up,is_start= relation['start1']<relation['end1'])
    vertex2 = chros[chr2].coordinate(relation['start2'],not is_up,is_start= relation['start2']<relation['end2'])
    vertex3 = chros[chr2].coordinate(relation['end2'],not is_up,is_start = relation['end2']<relation['start2'])
    vertex4 = chros[chr1].coordinate(relation['end1'],is_up,is_start = relation['end1']<relation['start1'])
    #ertex1,vertex2,vertex3,vertex4=chros[chr1].coordinate(relation['start1'],is_up,is_start=True),chros[chr2].coordinate(relation['start2'],not is_up,is_start=True),chros[chr2].coordinate(relation['end2'],not is_up),chros[chr1].coordinate(relation['end1'],is_up)
    svg+= '''<path d="M{} {} L{} {} L{} {} L{} {} Z" fill="{}" opacity="{}" />\n'''.format(vertex1[0],vertex1[1],vertex2[0],vertex2[1],vertex3[0],vertex3[1],vertex4[0],vertex4[1],relation['color'],relation['opacity'])
scale_x=args.svg_width*(1-args.svg_space/2)-L*scale;scale_y=args.svg_height*0.9
svg+='''<polyline points="{},{} {},{} {},{} {},{} " fill="none" stroke="black" stroke-width="3" />\n'''.format(scale_x,scale_y-5,scale_x,scale_y,scale_x+L*scale,scale_y,scale_x+L*scale,scale_y-5)
svg+='''<text x="{}" y="{}" fill="black"  font-size="{}">{}</text>\n'''.format(scale_x+L*scale/3,scale_y-10,args.svg_font_size,S)
if(args.highlight and os.path.exists(args.highlight)):
    h=open(args.highlight,'r')
    for line in h:
        if line.startswith('#'):continue
        line=line.strip()
        lines=re.split(r'\s+',line)
        if(len(lines)>=4):color=lines[3]
        else:color='red'
        c=lines[0];start=int(lines[1]);end=int(lines[2])+1
        svg+='''<rect x="{}" y="{}" width="{}" height="{}" fill="{}" />\n'''.format(chros[c].coordinate(start,is_up,is_start=True)[0],chros[c].top,chros[c].coordinate(end,is_up)[0]-chros[c].coordinate(start,is_up,is_start=True)[0],args.svg_thickness,color)
    h.close()


svg+='</svg>\n'
#=========================================================================

print(svg)
O=open(args.output+'.svg','w')
O.write(svg)
O.close()

os.system('convert {}.svg {}.png'.format(args.output,args.output))




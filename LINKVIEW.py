#!/usr/bin/env python
## _*_coding:utf-8_*_
from __future__ import division
import argparse
import re
from fractions import Fraction
import sys
import os
import interval


class ArgumentError(Exception):
    def __init__(self,ErrorType,ErrorValue):
        self.ErrorType = ErrorType
        self.ErrorValue = ErrorValue
    def __str__(self):
        if (self.ErrorType=="type"):
            return("[Error] Unkonwn input type '{}', only be 0 or 1 !".format(self.ErrorValue))
        if (self.ErrorType=="chr_len"):
            return("[Error] chr_len file '{}' not exists !".format(self.ErrorValue))
        if (self.ErrorType=="karyotype"):
            return("[Error] karyotype file '{}' not exists !".format(self.ErrorValue))
        if (self.ErrorType=="Parameter"):
            return("[Error] Parameter file '{}' not exists !".format(self.ErrorValue))
        if (self.ErrorType=="highlight"):
            return("Error] highlight file '{}' not exists !".format(self.ErrorValue))
        
class FormatError(Exception):
    def __init__(self,ErrorType,FileName,Line):
        self.ErrorType = ErrorType
        self.FileName = FileName
        self.Line = Line
    def __str__(self):
        if (self.ErrorType=="input"):
            return("\n\t[Error] Wrong format in karyotype file '{}':\n\t{}\n\tType parameter error or format error, please check carefully !".format(self.FileName,self.Line))
        if (self.ErrorType=="karyotype"):
            return("\n\t[Error] Wrong format in karyotype file '{}':\n\t{}\n\tThe correct format should be chr1:[start1:end1]...".format(self.FileName,self.Line))
        if (self.ErrorType=="karyotype_start_end"):
            return("\n\t[Error] Wrong format in karyotype file '{}':\n\t{}\n\tEnd must be larger than start".format(self.FileName,self.Line))

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

def main(args):
    
    global scale
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
            
        else: raise ArgumentError('type',args.type)
        if not opacity:opacity=args.opacity
        if not color:color=args.color
        try:
            start1=int(start1);end1=int(end1);start2=int(start2);end2=int(end2);
        except: raise FormatError('input',args.input,line)
        posofchro.setdefault(chr1,[]);posofchro.setdefault(chr2,[]);
        posofchro[chr1]+=[start1,end1];posofchro[chr2]+=[start2,end2];
        relations.append({'chr1':chr1,'start1':start1,'end1':end1,'chr2':chr2,'start2':start2,'end2':end2,'color':color,'opacity':opacity,'display':True})
    if(args.chr_len):
        if not os.path.exists(args.chr_len): raise ArgumentError('chr_len',args.chr_len)
        chr_len={}
        f=open(args.chr_len,'r')
        for line in f:
            line=line.strip()
            chr_name,length=re.match(r'^(\S+)\s+(\S+)\s*$',line).group(1,2)
            chr_len[chr_name]=int(length)
        f.close()


    chro_lst={}
    for c in posofchro:
        chro_lst[c]=chro(c,min(posofchro[c]),max(posofchro[c]),[True,True])
    order={}
    auto_start_end=[]
    if(args.karyotype):
        if not os.path.exists(args.karyotype): raise ArgumentError('karyotype',args.karyotype)
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
                if not name in chro_lst: sys.stderr.write('[Warning] {} has no alignment relationship with other chromosomes !\n'.format(name))
                if not len(arr) in (1,3): raise FormatError('karyotype',args.karyotype,line)
                if len(arr)==1:auto_start_end.append(name)
                if len(arr)==3:
                    if not (re.match(r'^\d+$',arr[1]) and re.match(r'^\d+$',arr[2])):
                        raise FormatError('karyotype',args.karyotype,line)
                    if not int(arr[1])<int(arr[2]):
                        raise FormatError('karyotype_start_end',args.karyotype,line)
                    chro_lst[name]=chro(name,int(arr[1]),int(arr[2]),[True,True])
                order[level].append(name)
            level+=1

        for i in range(len(relations)):
            relation=relations[i]
            chr1=relation['chr1'];chr2=relation['chr2']
            # chr1 在order中不存在 ，则这个relation 就不 display
            if not any([chr1 in x for x in order.values()]):relations[i]['display'] = False;continue
            if not any([chr2 in x for x in order.values()]):relations[i]['display'] = False;continue
            
            reverse = False
            if(relation['start1']-relation['end1'])*(relation['start2']-relation['end2'])<0: reverse=True;
            tmp_its1=interval.intersection([relation['start1'],relation['end1']] , [chro_lst[relation['chr1']].start,chro_lst[relation['chr1']].end])
            if not tmp_its1: relations[i]['display'] = False;continue;
            tmp_start1,tmp_end1 = sorted([relation['end1'],relation['start1']])
            len1 = tmp_end1-tmp_start1+1
            tmp_its1=[Fraction((tmp_its1[0]-tmp_start1),len1),Fraction((tmp_its1[1]-tmp_start1+1),len1)]
            tmp_its2=interval.intersection([relation['start2'],relation['end2']] , [chro_lst[relation['chr2']].start,chro_lst[relation['chr2']].end])
            if not tmp_its2: relations[i]['display']=False;continue;
            tmp_start2,tmp_end2 = sorted([relation['end2'],relation['start2']])
            len2 = tmp_end2-tmp_start2+1
            tmp_its2=[Fraction((tmp_its2[0]-tmp_start2),len2),Fraction((tmp_its2[1]-tmp_start2+1),len2)]
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
            chro_lst[c].start=min(posofchro[c])
            chro_lst[c]=chro(c,min(posofchro[c]),max(posofchro[c]),[True,True])

    else:
    #====Arrange the order of chro_lst on the image. If you specify a Karyotype file, this step is omitted.===
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
            nonlocal  chro_ordered
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
    svg_axis={}
    Gap_length={}

    max_len=0
    lens={}
    middle_spaces={}
    for level,chrs in order.items():
        num=len(chrs)
        length=0
        for c in chrs:length+=chro_lst[c].length
        svg_font_size[level]=args.svg_font_size
        show_pos_with_label[level]=args.show_pos_with_label
        svg_no_label[level]=args.svg_no_label
        svg_label_angle[level]=args.svg_label_angle
        svg_axis[level]=args.svg_axis
        Gap_length[level]=args.Gap_length
        lens[level]=length
    if(args.Parameter):
        if not os.path.exists(args.Parameter): raise ArgumentError('Parameter',args.Parameter)
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
            if m4:Gap_length[level]=float(m4.group(1))
            m5=re.search(r'svg_label_angle=(\d+(\.\d+)?)',line)
            if m5:svg_label_angle[level]=360-float(m5.group(1))
            m6 = re.match(r'svg_axis=(\d+(\.\d+)?)',line)
            if m6:svg_axis[level]=m6.group(1)
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
    L=max_len/5
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
    for level,chrs in order.items():
        left_start=(args.svg_width-lens[level]*scale)/2
        top+=space_vertical
        for c in chrs:
            if(args.chr_len and chro_lst[c].name in chr_len):
                chro_lst[c].is_full=[chro_lst[c].start==1,chro_lst[c].end==chr_len[chro_lst[c].name]]
            chro_lst[c].level=level
            chro_lst[c].top=top
            chro_lst[c].left=left_start
            chro_lst[c].right=left_start+chro_lst[c].length*scale
            if not chro_lst[c].is_full[0]: chro_lst[c].left+=line_length;chro_lst[c].right+=line_length
            left_start=chro_lst[c].right+middle_spaces[level]*scale
            if not chro_lst[c].is_full[1]: left_start+=line_length

    #==================creat svg=================

    def get_svg_content_Axis(c,position):
        axis_left = chro_lst[c].left
        axis_right = chro_lst[c].right
        axis_top = chro_lst[c].top
        axis_start=chro_lst[c].start
        axis_end = chro_lst[c].end
        
        for unit in [5,50,500,5000,50000,500000,5000000]:
            if args.svg_axis_density <= L/unit < 10*args.svg_axis_density:
                break
        unit2=unit
        unit1=int(unit2/5)
        unit3=unit2*2
        #print(unit1,unit2,unit3)
        def get_unit_start(unit,start):
            i = start
            while(1):
                if i%unit == 0: break
                i+=1
            unit_start=i
            return unit_start

        unit1_start = get_unit_start(unit1,axis_start)
        unit2_start = get_unit_start(unit2,axis_start)
        unit3_start = get_unit_start(unit3,axis_start)

        unit1_point=[]
        i=unit1_start
        while(i<=axis_end):
            unit1_point.append(i)
            i+=unit1

        unit2_point=[]
        i=unit2_start
        while(i<=axis_end):
            unit2_point.append(i)
            i+=unit2

        unit3_point=[]
        labels=[]
        i=unit3_start
        while(i<=axis_end):
            unit3_point.append(i)
            if unit3<=10000:
                labels.append('{:,.2f}'.format(i/1000)+' K')
            if unit3>10000:
                labels.append('{:,.2f}'.format(i/1000000)+' M')
            i+=unit3

        unit1_point = [x for x in unit1_point if not x in unit2_point]
        unit2_point = [x for x in unit2_point if not x in unit3_point]

        unit1_point = [(x-axis_start)*scale for x in unit1_point]
        unit2_point = [(x-axis_start)*scale for x in unit2_point]
        unit3_point = [(x-axis_start)*scale for x in unit3_point]

        axis_angle=25
        if position=='up':
            sign=-1
            axis_angle=360-axis_angle
        elif position=='down':
            sign=1
            axis_top+=args.svg_thickness
        svg_axis=''
        svg_axis+= '<path d="M {} {} L {} {} " stroke="black" fill="none"/>\n'.format(axis_left,axis_top,axis_right,axis_top)
        for point in unit1_point:
            svg_axis+='<path d="M {} {} L {} {} " stroke="black" fill="none"/>\n'.format(axis_left+point,axis_top,axis_left+point,axis_top+2*sign)
        for point in unit2_point:
            svg_axis+='<path d="M {} {} L {} {} " stroke="black" fill="none"/>\n'.format(axis_left+point,axis_top,axis_left+point,axis_top+5*sign)
        for i in range(len(unit3_point)):
            point=unit3_point[i]
            label=labels[i]
            svg_axis+='<path d="M {} {} L {} {} " stroke="black" fill="none"/>\n'.format(axis_left+point,axis_top,axis_left+point,axis_top+8*sign)
            svg_axis+='<text x="{}" y="{}" fill="black" transform="rotate({},{} {})" font-size="9">{}</text>\n'.format(axis_left+point,axis_top+12*sign,axis_angle,axis_left+point,axis_top+12*sign,label)
        return svg_axis

        
    max_level=len(order)-1

    svg_content_line = ''
    svg_content_chro = ''
    svg_content_label = ''
    svg_content_align = ''
    svg_content_scale = ''
    svg_content_highlight = ''
    svg_content_axis = ''


    for c in chro_lst:
        if not any([c in x for x in order.values()]):continue
        angle_text=''
        if not chro_lst[c].is_full[0]:
            svg_content_line+='''<line x1="{}" y1="{}" x2="{}" y2="{}" style="stroke:black;" class="line"/>\n'''.format(chro_lst[c].left-line_length,chro_lst[c].top+args.svg_thickness/2,chro_lst[c].left,chro_lst[c].top+args.svg_thickness/2)
        svg_content_chro+='''<rect x="{}" y="{}" width="{}" height="{}" fill="grey" class="chro"/>\n'''.format(chro_lst[c].left,chro_lst[c].top,chro_lst[c].right-chro_lst[c].left,args.svg_thickness)
        if not chro_lst[c].is_full[1]:
            svg_content_line+='''<line x1="{}" y1="{}" x2="{}" y2="{}" style="stroke:black;" class="line"/>\n'''.format(chro_lst[c].right,chro_lst[c].top+args.svg_thickness/2,chro_lst[c].right+line_length,chro_lst[c].top+args.svg_thickness/2)

        if(chro_lst[c].level==0):
            label_x=chro_lst[c].left;label_y=chro_lst[c].top-args.svg_thickness/2
            if svg_axis[chro_lst[c].level]:
                svg_content_axis += get_svg_content_Axis(c,'up')
                label_y -= 22
            
        elif(chro_lst[c].level==max_level):
            label_x=chro_lst[c].left;label_y=chro_lst[c].top+args.svg_thickness*2;
            if svg_axis[chro_lst[c].level]:
                svg_content_axis += get_svg_content_Axis(c,'down')
                label_y += 22
        else:
            label_x=chro_lst[c].right+svg_font_size[chro_lst[c].level]/5;
            label_y=chro_lst[c].top+args.svg_thickness
            if not chro_lst[c].is_full[1]: label_x+=line_length
        if not int(svg_no_label[chro_lst[c].level]):
            label=chro_lst[c].name
            if int(show_pos_with_label[chro_lst[c].level]):label="{} {:,} ~ {:,} bp".format(label,chro_lst[c].start,chro_lst[c].end);
            angle_text='transform="rotate({},{} {})"'.format(svg_label_angle[chro_lst[c].level],label_x,label_y)
            svg_content_label+='''<text x="{}" y="{}" fill="black" font-size="{}" {} class="label">{}</text>\n'''.format(label_x,label_y,svg_font_size[chro_lst[c].level],angle_text,label)
    for relation in relations:
        chr1=relation['chr1'];chr2=relation['chr2']
        is_up=chro_lst[chr1].level>chro_lst[chr2].level
        vertex1 = chro_lst[chr1].coordinate(relation['start1'],is_up,is_start= relation['start1']<relation['end1'])
        vertex2 = chro_lst[chr2].coordinate(relation['start2'],not is_up,is_start= relation['start2']<relation['end2'])
        vertex3 = chro_lst[chr2].coordinate(relation['end2'],not is_up,is_start = relation['end2']<relation['start2'])
        vertex4 = chro_lst[chr1].coordinate(relation['end1'],is_up,is_start = relation['end1']<relation['start1'])
        
        svg_content_align+= '''<path d="M{} {} L{} {} L{} {} L{} {} Z" fill="{}" opacity="{}" class="align"/>\n'''.format(vertex1[0],vertex1[1],vertex2[0],vertex2[1],vertex3[0],vertex3[1],vertex4[0],vertex4[1],relation['color'],relation['opacity'])
    scale_x=args.svg_width*(1-args.svg_space/2)-L*scale;scale_y=args.svg_height*0.9
    svg_content_scale+='''<polyline points="{},{} {},{} {},{} {},{} " fill="none" stroke="black" class="scale"/>\n'''.format(scale_x,scale_y-5,scale_x,scale_y,scale_x+L*scale,scale_y,scale_x+L*scale,scale_y-5)
    svg_content_scale+='''<text x="{}" y="{}" fill="black"  font-size="{}" class="scale-text">{}</text>\n'''.format(scale_x+L*scale/3,scale_y-10,args.svg_font_size,S)
    if(args.highlight):
        if not os.path.exists(args.highlight): raise ArgumentError('highlight',args.highlight)
        h=open(args.highlight,'r')
        for line in h:
            if line.startswith('#'):continue
            line=line.strip()
            lines=re.split(r'\s+',line)
            if(len(lines)>=4):color=lines[3]
            else:color='red'
            c=lines[0];start=int(lines[1]);end=int(lines[2])+1
            svg_content_highlight+='''<rect x="{}" y="{}" width="{}" height="{}" fill="{}" />\n'''.format(chro_lst[c].coordinate(start,is_up,is_start=True)[0],chro_lst[c].top,chro_lst[c].coordinate(end,is_up)[0]-chro_lst[c].coordinate(start,is_up,is_start=True)[0],args.svg_thickness,color)
        h.close()


    svg_content = '<svg width="{}" height="{}" xmlns="http://www.w3.org/2000/svg" version="1.1">\n'.format(args.svg_width,args.svg_height)
    svg_content += (svg_content_line + svg_content_chro + svg_content_label + svg_content_align + svg_content_scale + svg_content_highlight+svg_content_axis)
    svg_content += '</svg>\n'
    #=========================================================================

    print(svg_content)
    O=open(args.output+'.svg','w')
    O.write(svg_content)
    O.close()

    os.system('convert {}.svg {}.png'.format(args.output,args.output))   

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input",help='The input file which contains location relationships. Can be a text file, format like: "chr1 start1 end1 chr2 start2 end2 color:opacity". It can also be a blast output file.')
    parser.add_argument('-t','--type',default=0,type=int,help='Type of input file, 0 for link.txt format like "chr1 start1 end1 chr2 start2 end2 color:opacity", 1 for output file of blast format as "query subject identity alignment_length mismatches gap_opens q.start q.end s.start s.end evalue bit score"\tdefault=0')
    parser.add_argument('-hl','--highlight',help='highlight.txt, Highlight section, format like chr start end [color]')
    parser.add_argument('-k','--karyotype',help="Karyotype.txt, If you don't specify this file, it will automatically sort and display.format like chr1:[start1:end1] chr2...")
    parser.add_argument('--opacity',default=0.5,type=float,help='opacity, default=0.5')
    parser.add_argument('--color',default="green",type=str,help='color, default=green')
    parser.add_argument('--svg_height',default=800,type=int,help='height of svg, default=800')
    parser.add_argument('--svg_width',default=1200,type=int,help='width of svg, default=1200')
    parser.add_argument('--svg_thickness',default=15,type=int,help='thickness of chromosome, default=15')
    parser.add_argument('-n','--svg_no_label',action="store_true",help='Do not show labels')
    parser.add_argument('--svg_font_size',default=18,type=int,help='font size of the label, default=18')
    parser.add_argument('--svg_label_angle',default=0,type=int,help='label rotation angle, default=0')
    parser.add_argument('--svg_axis',action="store_true",help='Display the axis (only for the chromosomes at the top or bottom)')
    parser.add_argument('--svg_axis_density',default=2,type=int,help='It is useful for tuning the scale density. The higher the value, the denser the scale.')
    parser.add_argument('--svg_space',default=0.2,type=float,help='The proportion of white space left and right, default=0.2')
    parser.add_argument('-s','--show_pos_with_label',action="store_true",help='Display location information on the label')
    parser.add_argument('--scale',help='example:5k default=auto')
    parser.add_argument('-o','--output',default="linkview_output",type=str,help='output file prefix, default=linkview_output')
    parser.add_argument('--min_identity',default=95,type=float,help='if your input file is a output file of blast, the min_identy to filter the input file, default=95')
    parser.add_argument('--min_alignment_length',default=200,type=int,help='if your input file is a output file of blast, the min_alignment_length to filter the input file, default=200')
    parser.add_argument('--max_evalue',default=1e-5,type=float,help='if your input file is a output file of blast, the max_evalue to filter the input file, default=1e-5')
    parser.add_argument('--min_bit_score',default=5000,type=float,help='if your input file is a output file of blast, the min_bit_score to filter the input file, default=5000')
    parser.add_argument('--Gap_length',default=0.2,type=float,help='Length of gap between two segments per line,if > 1,It represents Physical length, if<1,It represents total_length_of_this_line * this. default=0.2')
    parser.add_argument('--chr_len',help='chr_len.txt, format as: chr1 150000')
    parser.add_argument('-P','--Parameter',help='Specify the parameters for each row separately in a file, if needed. E.g svg_font_size=20 show_pos_with_label=0 svg_no_label=1 Gap_length=500')
    args=parser.parse_args()
    if args.svg_label_angle:args.svg_label_angle=360-float(args.svg_label_angle)

    main(args)


# LINKVIEW
LINKVIEW 是一个将序列比对结果（或其它任何适合的数据）进行可视化作图的工具。
本工具设计灵感来源于circos软件，因为在日常工作中需要对blast比对结果进行可视化显示，但是找不到合适的软件(circos仅能作环状的图)，所以开发了这款工具。
### python版本
在 2.7.9 和 3.6.2 测试通过
### USAGE
LINKVIEW 的USAGE如下：
```
usage: LINKVIEW.py [-h] [-t TYPE] [-hl HIGHLIGHT] [-k KARYOTYPE]
                   [--opacity OPACITY] [--color COLOR]
                   [--svg_height SVG_HEIGHT] [--svg_width SVG_WIDTH]
                   [--svg_thickness SVG_THICKNESS] [-n]
                   [--svg_font_size SVG_FONT_SIZE]
                   [--svg_label_angle SVG_LABEL_ANGLE] [--svg_axis]
                   [--svg_axis_density SVG_AXIS_DENSITY]
                   [--svg_space SVG_SPACE] [-s] [--scale SCALE] [-o OUTPUT]
                   [--min_identity MIN_IDENTITY]
                   [--min_alignment_length MIN_ALIGNMENT_LENGTH]
                   [--max_evalue MAX_EVALUE] [--min_bit_score MIN_BIT_SCORE]
                   [--Gap_length GAP_LENGTH] [--chr_len CHR_LEN]
                   [-P PARAMETER]
                   input

positional arguments:
  input                 The input file which contains location relationships.
                        Can be a text file, format like: "chr1 start1 end1
                        chr2 start2 end2 color:opacity". It can also be a
                        blast output file.

optional arguments:
  -h, --help            show this help message and exit
  -t TYPE, --type TYPE  Type of input file, 0 for link.txt format like "chr1
                        start1 end1 chr2 start2 end2 color:opacity", 1 for
                        output file of blast format as "query subject identity
                        alignment_length mismatches gap_opens q.start q.end
                        s.start s.end evalue bit score" default=0
  -hl HIGHLIGHT, --highlight HIGHLIGHT
                        highlight.txt, Highlight section, format like chr
                        start end [color]
  -k KARYOTYPE, --karyotype KARYOTYPE
                        Karyotype.txt, If you don't specify this file, it will
                        automatically sort and display.format like
                        chr1[:start1:end1] chr2...
  --opacity OPACITY     opacity, default=0.5
  --color COLOR         color, default=green
  --svg_height SVG_HEIGHT
                        height of svg, default=800
  --svg_width SVG_WIDTH
                        width of svg, default=1200
  --svg_thickness SVG_THICKNESS
                        thickness of chromosome, default=15
  -n, --svg_no_label    Do not show labels
  --svg_font_size SVG_FONT_SIZE
                        font size of the label, default=18
  --svg_label_angle SVG_LABEL_ANGLE
                        label rotation angle, default=0
  --svg_axis            Display the axis (only for the chromosomes at the top
                        or bottom)
  --svg_axis_density SVG_AXIS_DENSITY
                        It is useful for tuning the scale density. The higher
                        the value, the denser the scale.
  --svg_space SVG_SPACE
                        The proportion of white space left and right,
                        default=0.2
  -s, --show_pos_with_label
                        Display location information on the label
  --scale SCALE         example:5k default=auto
  -o OUTPUT, --output OUTPUT
                        output file prefix, default=linkview_output
  --min_identity MIN_IDENTITY
                        if your input file is a output file of blast, the
                        min_identy to filter the input file, default=95
  --min_alignment_length MIN_ALIGNMENT_LENGTH
                        if your input file is a output file of blast, the
                        min_alignment_length to filter the input file,
                        default=200
  --max_evalue MAX_EVALUE
                        if your input file is a output file of blast, the
                        max_evalue to filter the input file, default=1e-5
  --min_bit_score MIN_BIT_SCORE
                        if your input file is a output file of blast, the
                        min_bit_score to filter the input file, default=5000
  --Gap_length GAP_LENGTH
                        Length of gap between two segments per line,if > 1,It
                        represents Physical length, if<1,It represents
                        total_length_of_this_line * this. default=0.2
  --chr_len CHR_LEN     chr_len.txt, format as: chr1 150000
  -P PARAMETER, --Parameter PARAMETER
                        Specify the parameters for each row separately in a
                        file, if needed. E.g svg_font_size=20
                        show_pos_with_label=0 svg_no_label=1 Gap_length=500
```
### 使用方法说明

1. 输入文件 <br>
    在使用LINKVIEW作图前，需要先准备一个如下格式的文件：
    
    ```chr1 start1 end1 chr2 start2 end2 color:opacity```

    chr1和chr2表示序列(染色体)名称；
    start1、end1、start2、end2为位点位置(bp，整数)，表示chr1的start1到end1比对上chr2的start2到end2；
    color:opacity 表示该比对块绘制的颜色和透明度，可以省略，缺省值为green:0.5。
    该文件每行为一个alignment。
    可以将任何比对软件的结果(或其他任何表示对应关系的数据)转换成如上格式，即可通过LINKVIEW作图显示。


    输入文件也可以为blastn的输出结果(tabular格式): 
    ```query id, subject id, % identity, alignment length, mismatches，gap opens, q. start, q. end, s. start, s. end, evalue, bit score```
    此时需要指定--type 1

    当准备好输入文件后，即可运行LINKVIEW：
    ./LINKVIEW.py [-t TYPE] input
    LINKVIEW 会自动分配染色体在图中的位置。
    <br>
2. KARYOTYPE文件 <br>
    若不指定KARYOTYPE文件，LLINKVIEW会绘制输入文件中的所有染色体，并自动分配它们在图中的位置。
    通过-k参数指定一个KARYOTYPE文件，可以指定需要绘制的染色体及其在图中的位置和所需绘制的区间。
    KARYOTYPE 格式下：
    ```
    chr1[:start1:end1] chr2[:start2:end2]
    chr3[:start3:end3]
    ```
    每一行对应图中的每一横排，上面内容的含义是将chr1和chr2绘制在第一横排，且chr1在左，chr2在右，chr3绘制在第二横排；
    start和end指定需要绘制的该染色体的区间，可以省略，若省略，则LINKVIEW会根据输入文件计算出包含所有alignment的最小区间。
    <br>
3. HIGHLIGHT文件 <br>
    LINKVIEW 可以将染色体上部分区段高亮显示。
    通过 -hl参数指定一个HIGHLIGHT文件，格式如下：
    ```chr start end [color]```
    该文件每一行表示一个高亮区块，上面内容表示chr的start到end区间显示为高亮，
    可以在最后指定颜色，否则为默认颜色为red
    <br>
4. CHR_LEN文件 <br>
   通过--chr_len参数指定CHR_LEN文件，格式如下：
    ```
    chr1 len1
    chr2 len2
    ```
    该文件指定每条染色体的长度，如果指定该文件，KARYOTYPE文件中指定的染色体区间不完整时，不完整的部分将以一条短横线代替
    <br>
5. 其它参数 <br>
   
  <table>
  <tr>
		<th>-o, --output</th>
		<td>输出文件的前缀</td>
	</tr>
	<tr>
    <th>-n, --svg_no_label</th>
		<td>不显示标签</td>
	</tr>
	<tr>
    <th>--svg_font_size</th>
    <td>标签字体大小</td>
  </tr>
	<tr>
    <th>--svg_label_angle</th>
    <td>标签旋转角度</td>
  </tr>
	<tr>
    <th>--svg_axis</th>
    <td>显示刻度(仅对第一横排和最后一横排有效)</td>
  </tr>
	<tr>
    <th>--svg_axis_density</th>
    <td>大于0的值，对于控制刻度的稠密程度有一定作用</td>
  </tr>
  <tr>
    <th>--show_pos_with_label</th>
    <td>显示标签的同时，也显示染色体的区间位置信息</td>
  </tr>
  <tr>
    <th>--scale</th>
    <td>在右下角绘制的比例尺的大小，默认为自动大小</td>
  </tr>
  <tr>
    <th>--min_identity</th>
    <td rowspan="43">当输入文件类型为blastn结果时，用于过滤alignment</td>
  </tr>
  <tr>
    <th>--min_alignment_length</th>
    
  </tr>
  <tr>
    <th>--max_evalue</th>
    
  </tr>
  <tr>
    <th>--min_bit_score</th>
    
  </tr>
</table>

值得注意的是，可以通过-P指定PARAMETER文件，为每一横排指定参数,PARAMETER文件的每一行对应作图的每一横排，格式为parameter=value，例如：
```
svg_axis=1 svg_label_angle=30
svg_font_size=20
```
上面内容的含义是，作图的第一横排显示刻度，标签旋转30°，第二横排字体大小为20px

<br>
<br>
<hr>

联系邮箱：<br>
<a href="mailto:shunlintianxia@sina.com">shunlintianxia@sina.com</a><br>
<a href="mailto:397441459@qq.com">397441459@qq.com</a>


         


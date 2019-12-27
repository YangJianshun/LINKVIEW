# LINKVIEW
LINKVIEW 是一个将序列比对结果（或其它任何适合的数据）进行可视化作图的工具。
本工具设计灵感来源于circos软件，因为在日常工作中需要对blast比对结果进行可视化显示，但是找不到合适的软件(circos仅能作环状的图)，所以开发了这款工具。
使用LINKVIEW可以绘制出整体或局部的比对情况，支持自定义高亮、导入gff文件以绘制基因结构、有多种风格供选择。请看效果图：
![染色体比对图](imgs/example2.png)
![局部比对图](imgs/example3.png)

### 项目依赖

#### python
在 python 2.7.9 和 3.6.2 测试通过
#### inkscape
svg转为png时需要调用inkscape，在 inkscape 0.91 测试通过

### 使用方法说明

1. 输入文件 <br>
    在使用LINKVIEW作图前，需要先准备一个alignment文件，可以自定义，也可以直接使用blastn或MUMmer的比对结果。
    自定义的输入文件格式如下：
    ```chr1 start1 end1 chr2 start2 end2 [color:opacity]```
    该文件每行为一个alignment。
    chr1和chr2表示序列(染色体)名称；
    start1、end1、start2、end2为位点位置(bp，整数)，表示chr1的start1到end1比对上chr2的start2到end2；
    color:opacity 表示该比对块绘制的颜色和透明度，可以省略，缺省值按照指定的风格有所不同。
    <br>

    可以之间将比对软件的输出作为LINKVIEW的输入文件，目前支持blastn(tabular格式，见example3)和MUMmer(show-coords生成的文件，见example1)

    <br>

    当准备好输入文件后，即可运行LINKVIEW：
    ```./LINKVIEW.py [-t TYPE] input```

    <br>
2. KARYOTYPE文件 <br>
    若不指定KARYOTYPE文件，LINKVIEW会绘制输入文件中的所有染色体，并自动分配它们在图中的位置。
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
    颜色可以省略，默认颜色为red
    <br>
4. CHR_LEN文件 <br>
   通过--chr_len参数指定CHR_LEN文件，格式如下：
    ```
    chr1 len1
    chr2 len2
    ```
    该文件指定每条染色体的长度，如果指定该文件，KARYOTYPE文件中指定的染色体区间不完整展示时，没展示的部分将以一条短横线代替
    <br>
5. 其它参数 <br>
   
  <table>
  <tr>
		<th>-o, --output</th>
		<td>输出文件的前缀</td>
	</tr>
	<tr>
    <th>-n, --no_label</th>
		<td>不显示标签</td>
	</tr>
	<tr>
    <th>--label_font_size</th>
    <td>标签字体大小</td>
  </tr>
	<tr>
    <th>--label_angle</th>
    <td>标签旋转角度</td>
  </tr>
	<tr>
    <th>--chro_axis</th>
    <td>显示刻度(仅对第一横排和最后一横排有效)</td>
  </tr>
	<tr>
    <th>--chro_axis_density</th>
    <td>大于0的值，对于控制刻度的稠密程度有一定作用</td>
  </tr>
  <tr>
    <th>--show_pos_with_label</th>
    <td>显示标签的同时，也显示染色体的区间位置信息</td>
  </tr>
  <tr>
    <th>--scale</th>
    <td>在右下角绘制的比例尺线段图例的大小，默认为自动大小</td>
  </tr>
  <tr>
    <th>--gff</th>
    <td>指定gff文件以绘制基因结构</td>
  </tr>
  <tr>
    <th>--bezier</th>
    <td>比对块绘制成贝塞尔曲线风格</td>
  </tr>
  <tr>
    <th>--style</th>
    <td>改变绘图风格，目前有两种样式分别是classic和simple，前者颜色较深，后者较浅</td>
  </tr>
  <tr>
    <th>--min_identity</th>
    <td rowspan="43">用于过滤alignment</td>
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

值得注意的是，可以通过-P指定PARAMETER文件，为每一横排指定参数，PARAMETER文件的每一行对应作图的每一横排，格式为parameter=value，例如：
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


         


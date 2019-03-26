#raw=[[1,7],[2,5],[6,9],[13,17],[12,18],[101,102]]
#__all__ = ('forward','reverse','merge','relation')
__all__ = ('merge','relation')
def forward(intervals):
    return [sorted(x) for x in intervals]
def reverse(intervals):
    return [sorted(x,reverse = True) for x in intervals]


def merge(intervals):
    intervals=forward(intervals)
    intervals.sort(key=lambda x: x[0])
    
    merged = []
    for interval in intervals:
        if not merged or merged[-1][1] < interval[0]: 
            merged.append(interval)
        else:
             merged[-1][1] = max(merged[-1][1], interval[1])
    return merged

def relation(interval1,interval2):
    min1,max1=sorted(interval1)
    min2,max2=sorted(interval2)
    if(min1==min2 and max1==max2):return 0
    if(max1<min2 or max2<min1):return 1
    if(min1<min2<=max1<max2 or min2<min1<=max2<max1):return 2
    if(min1<=min2<=max2<=max1 or min2<=min1<=max1<=max2):return 3
    #0 xiang deng
    #1 xiang li
    #2 xiang jiao
    #3 bao han

def is_in(interval1,interval2):
    min1,max1=sorted(interval1)
    min2,max2=sorted(interval2)
    if(min2<=min1<=max1<=max2):return True
    else: return False
   
def complement(wide,intervals):
    s,e=wide
    front,back=s,s
    result=[]
    intervals=forward(intervals)
    intervals.sort(key=lambda x: x[0])
    for (s1,e1) in intervals:
        front=s1
        if front>back:result.append([back,front-1])
        if e1+1>back: back=e1+1
    if back<=e: result.append([back,e])
    return(result)
def intersection(interval1,interval2):
    nums=sorted(interval1+interval2)
    if relation(interval1,interval2) != 1:
        return [nums[1],nums[2]]
    else:return False
    

#pos=[[2,3],[7,-4]]
#pos=[[-10,-7],[-3,-1],[8,9],[45,108]]
#total=[0,100]
#total=[0, 11]
#pos=[[-7,8],[0, 2], [6, 7],[10,100]]
#print(complement(total,pos))

#after=merge(raw)
#print(after)

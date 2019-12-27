__all__ = ('merge','relation','is_in','complement','intersection')

def merge(intervals):
    '''
    @msg: Merge multiple intervals
    @param intervals {list} A two-dimensional array, each item representing an interval
    @return: {list}  Return the merged interval list
    '''

    intervals = [sorted(x) for x in intervals]
    intervals.sort(key=lambda x: x[0])
    merged = []
    for interval in intervals:
        if not merged or merged[-1][1] < interval[0]:
            merged.append(interval)
        else:
             merged[-1][1] = max(merged[-1][1], interval[1])
    return merged

def relation(interval1,interval2):
    '''
    @msg: Judging the relationship between two intervals
    @param interval1 {list} The first interval
    @param interval2 {list} The second interval
    @return: {int}  Return the relationship between two intervals, 0: the two intervals are equal, 1: the two intervals are separated, 2: the two intervals intersect, 3: one interval is contained by the other interval
    '''

    min1,max1=sorted(interval1)
    min2,max2=sorted(interval2)
    if(min1==min2 and max1==max2):return 0
    if(max1<min2 or max2<min1):return 1
    if(min1<min2<=max1<max2 or min2<min1<=max2<max1):return 2
    if(min1<=min2<=max2<=max1 or min2<=min1<=max1<=max2):return 3

def is_in(interval1,interval2):
    '''
    @msg: Determine if an interval is contained by another interval
    @param interval1 {list} The first interval
    @param interval2 {list} The second interval
    @return: {bool}  interval1 contained by interval2 return true, otherwise return false
    '''

    min1,max1=sorted(interval1)
    min2,max2=sorted(interval2)
    if(min2<=min1<=max1<=max2):return True
    else: return False

def complement(wide,intervals):
    '''
    @msg: Find the complement of multiple intervals over a large interval
    @param wide {list} The large interval
    @param intervals {list} List of multiple intervals
    @return: {list}  Return a list of complement intervals
    '''

    start,end=wide
    front,back=start,start
    result=[]
    intervals = [sorted(x) for x in intervals]
    intervals.sort(key=lambda x: x[0])
    for (start1,end1) in intervals:
        if start1 > end: break
        front = start1
        if front > back: result.append([back,front-1])
        if end1+1 > back: back=end1+1
    if back<=end: result.append([back,end])
    return(result)

def intersection(interval1,interval2):
    '''
    @msg: Find the intersection of two intervals
    @param interval1 {list} The first interval
    @param interval2 {list} The second interval
    @return: {list}  If there is an intersection between two intervals, return the intersection, otherwise return an empty list
    '''

    nums=sorted(interval1+interval2)
    if relation(interval1,interval2) != 1:
        return [nums[1],nums[2]]
    else:return []



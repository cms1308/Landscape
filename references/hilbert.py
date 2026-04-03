import json
import sys
import requests
import subprocess
import multiprocessing as mp
import ast
import time
from itertools import product
from itertools import repeat
from multiprocessing import freeze_support
import os

direc = os.getcwd()
# def tilde(str):
#     result='Overtilde[%s]'%(str)
#     return result
#
# def adj():
#     return '\[CaptialPhy]'




def group(rank):
    result = 'C'+str(rank)
    return result


def matter(rank, num):
    if rank >= 2:
        xsinglet = [0]*rank
        msinglet = [0]*rank
        fund = [0] * rank
        fund[0] = 1
        antifund = [0] * rank
        antifund[-1] = 1
        adj = [0] * rank
        adj[0] = 1
        adj[-1] = 1
        sym = [0] * rank
        sym[0] = 2
        symbar = [0] * rank
        symbar[-1] = 2
        anti = [0] * rank
        anti[1] = 1
        antibar = [0] * rank
        antibar[-2] = 1

    else:
        xsinglet = [0]
        msinglet = [0]
        fund = [1]
        antifund = [1]
        adj = [2]
        sym = [2]
        symbar = [2]
        anti = [0]
        antibar = [0]
    list = [xsinglet, msinglet, fund, antifund, adj, sym, symbar, anti, antibar]
    result = []
    for i in range(len(num)):
        for j in range(num[i]):
            result.append(list[i])
    return str(result).replace(" ", "")


def matters2(rank, num):
    if rank >= 2:
        xsinglet = [0]*rank
        msinglet = [0]*rank
        fund = [0] * rank
        fund[0] = 1
        antifund = [0] * rank
        antifund[-1] = 1
        adj = [0] * rank
        adj[0] = 1
        adj[-1] = 1
        sym = [0] * rank
        sym[0] = 2
        symbar = [0] * rank
        symbar[-1] = 2
        anti = [0] * rank
        anti[1] = 1
        antibar = [0] * rank
        antibar[-2] = 1

    else:
        xsinglet = [0]
        msinglet = [0]
        fund = [1]
        antifund = [1]
        adj = [2]
        sym = [2]
        symbar = [2]
        anti = [0]
        antibar = [0]
    list = [xsinglet, msinglet, fund, antifund, adj, sym, symbar, anti, antibar]
    result=[]
    for i in range(len(num)):
        if num[i] != 0:
            result.append(list[i])
    return str(result).replace(" ", "")


def num(nums):
    result = [1]*sum(nums)
    return str(result).replace(" ", "")


def num2(nums):
    result=[]
    for i in range(len(nums)):
        if nums[i]!=0:
            result.append(nums[i])
    return str(result).replace(" ","")

############## without multiprocessing ###############
# part=str(list(product(range(1,2),range(0,1),range(0,1),range(0,1),range(0,1)))).replace("(","[").replace(")","]")
#
# lcode = """
#     maxnodes 9999999
#     maxobjects 999999
#     group=%s
#     matter=%s
#     num=%s
#     part=%s
#     comb(int n,m)={
#     class_ord([n+1])/(class_ord([n-m+1])*class_ord([m+1]))
#     }
#     single=poly_one(n_cols(matter))
#     sol=null(0,n_cols(part)+1)
#         for j=1 to n_rows(part) do
#             term=single;
#             for k=1 to size(part[j]) do
#                 if part[j,k]==0 then term=tensor(term,single,group);
#                 else frob=from_part(partitions(part[j,k]));
#                     ffrob=null(0,n_cols(frob)+1);
#                     for l=1 to n_rows(frob) do
#                         if all_one(part[j,k])*(frob[l]+(l/(n_rows(frob))))<=num[k] then ffrob=ffrob+(frob[l]+(l/(n_rows(frob)))); fi;
#                     od;
#                 zero=0*single;
#                 for m=1 to n_rows(ffrob) do
#                     rep=single;
#                     count=num[k];
#                     for n=1 to n_cols(ffrob) do
#                         rep=comb(count,ffrob[m,n])*tensor(rep,p_tensor(ffrob[m,n], sym_tensor(n,matter[k],group),group),group);
#                         count=count-ffrob[m,n];
#                     od;
#                     zero=zero+rep;
#                 od;
#                 term=tensor(term,zero,group);
#                 fi;
#             od;
#             if coef(term,1)!=0 && term[1]/coef(term,1)==single then sol=sol+(part[j]+coef(term,1)); fi;
#         od;
#     print(sol);
#
#     """ % (groups, matters, nums, part)
# res = subprocess.run(['lie'], input=lcode, capture_output=True, encoding='UTF-8').stdout[53:].strip()


def PE(count, matters, parts, nums, rank):
    lcode="""
    maxnodes 9999999
    maxobjects 9999999
    group=%s
    matter=%s
    num=%s
    part=%s
    comb(int n,m)={
    class_ord([n+1])/(class_ord([n-m+1])*class_ord([m+1]))
    }
    single=poly_one(n_cols(matter))
    sol=null(0,n_cols(part)+1)
        for j=1 to n_rows(part) do
            term=single;
            for k=1 to size(part[j]) do
                if part[j,k]==0 then term=tensor(term,single,group);
                else frob=from_part(partitions(part[j,k]));
                    ffrob=null(0,n_cols(frob)+1);
                    for l=1 to n_rows(frob) do
                        if all_one(part[j,k])*(frob[l]+(l/(n_rows(frob))))<=num[k] then ffrob=ffrob+(frob[l]+(l/(n_rows(frob)))); fi;
                    od;
                zero=0*single;
                for m=1 to n_rows(ffrob) do
                    rep=single;
                    count=num[k];
                    for n=1 to n_cols(ffrob) do
                        rep=comb(count,ffrob[m,n])*tensor(rep,p_tensor(ffrob[m,n], sym_tensor(n,matter[k],group),group),group);
                        count=count-ffrob[m,n];
                    od;
                    zero=zero+rep;
                od;
                term=tensor(term,zero,group);
                fi;
            od;
            if coef(term,1)!=0 && term[1]/coef(term,1)==single then sol=sol+(part[j]+coef(term,1)); fi;
        od;
    print(sol);
    """ % (group(rank), matters, nums, '['+str(parts[count])+']')
    res = subprocess.run(['lie'], shell=True, input=lcode, capture_output=True, encoding='UTF-8').stdout[53:].replace('\n', '').replace(" ", "")
    if 'null' in res:
        pass
    else:
        return (res)


def PL(strs):
    mcode = """
    count = %s; count2 = {}; count3 = {};
    For[i = 1, i <= Length[count], i++, count4 = 1;
        For[j = 1, j <= count[[i]], j++, count2 = Append[count2, i];
            count3 = Append[count3, count4]; count4++;]];
    seriesvar = Subscript[z, count2[[#]], count3[[#]]] & /@ Range[Length[count2]];
    f = Total[#[[Length[#]]] Product[Subscript[z, count2[[i]], count3[[i]]]^#[[i]], {i, 1,Length[#] - 1}] & /@ % s];
    Print[f];
    PL[x_, k_] := Sum[(MoebiusMu[i] Log[x])/i /.Subscript[z,w_,y_] :>  Subscript[z,w, y]^i, {i, 1, k}];
    pl = PL[f, %d];
    For[i = 1, i <= Length[%s[[1]]]-1, i++,pl = Series[pl, {seriesvar[[i]], 0, %d}] // Normal // Expand;];
    Print[pl];
    """ % (counts, strs, order, strs, order)
    proc = subprocess.Popen(['wolframscript', '-code', mcode], stdout=subprocess.PIPE)
    (out, err) = proc.communicate()
    result = out.decode('ascii').replace("Null", "")
    # replace("Subscript","").replace(", ","_")
    return result


def PL2(strs):
    mcode = """
    count = %s; count2 = {}; count3 = {};count4 = 1;
    For[i = 1, i <= Length[count], i++,
        If[count[[i]]!=0, count2=Append[count2, count4]; count4++;]; 
    ];

    seriesvar = Subscript[z, count2[[#]]] & /@ Range[Length[count2]];
    f = Total[#[[Length[#]]] Product[Subscript[z, count2[[i]]]^#[[i]], {i, 1,Length[#] - 1}] & /@ % s];
    Print[f];
    PL[x_, k_] := Sum[(MoebiusMu[i] Log[x])/i /.Subscript[z,w_] :>  Subscript[z,w]^i, {i, 1, k}];
    pl = PL[f, %d];
    For[i = 1, i <= Length[%s[[1]]]-1, i++,pl = Series[pl, {seriesvar[[i]], 0, %d}] // Normal // Expand;];
    Print[pl];
    """ % (counts, strs, order, strs, order)
    proc = subprocess.Popen(['wolframscript', '-code', mcode], stdout=subprocess.PIPE)
    (out, err) = proc.communicate()
    result = out.decode('ascii').replace("Null", "")
    # replace("Subscript","").replace(", ","_")
    return result


def Hilbert(strs):
    mcode = """
        count = %s; count2 = {}; count3 = {};
        For[i = 1, i <= Length[count], i++, count4 = 1;
            For[j = 1, j <= count[[i]], j++, count2 = Append[count2, i];
                count3 = Append[count3, count4]; count4++;]];
        seriesvar = Subscript[z, count2[[#]], count3[[#]]] & /@ Range[Length[count2]];
        f = #[[Length[#]]] Product[Subscript[z, count2[[i]], count3[[i]]]^#[[i]], {i, 1,Length[#] - 1}] & /@ % s
        Print[f];ng1
        """ % (counts, strs)
    proc = subprocess.Popen(['wolframscript', '-code', mcode], stdout=subprocess.PIPE)
    (out, err) = proc.communite()
    result = out.decode('ascii').replace("Null", "")  # .replace("Subscript","").replace(", ","_")
    return result


def run(matters, parts, nums, rank):
    with mp.Pool() as pool:
        # ans = pool.map(PE, range(len(parts)))
        ans = pool.starmap(func=PE, iterable=zip(range(len(parts)), repeat(matters), repeat(parts), repeat(nums), repeat(rank)))
    list2 = list(filter(None, ans))
    pe=str(list2).replace("[[", "{").replace("]]", "}").replace("[", "{").replace("]", "}").replace("'", "").replace(" ", "")
    return PL(pe)


def run2(matters, parts, nums, rank):
    with mp.Pool() as pool:
        # ans = pool.map(PE, range(len(parts)))
        ans = pool.starmap(func=PE, iterable=zip(range(len(parts)), repeat(matters), repeat(parts), repeat(nums), repeat(rank)))
    list2 = list(filter(None, ans))
    pe=str(list2).replace("[[", "{").replace("]]", "}").replace("[", "{").replace("]", "}").replace("'", "").replace(" ", "")
    return PL2(pe)


if __name__ == '__main__':
    start = time.time()
    rank = 4
    count = [0, 0, 0, 0, 0, 1, 0, 1, 0] #[X,M,q,qb,adj,s,sb,a,ab]
    counts=str(count).replace("[", "{").replace("]", "}")
    order = 6
    nums = num2(count)
    direc = "/Users/cms1308/Library/CloudStorage/OneDrive-Personal/Graduate/shared folder/Mathematica/Hilbert series"
    
    matters = matters2(rank,count)
    parts = ast.literal_eval(
        str(list(product(range(0, order + 1), repeat=len(ast.literal_eval(nums))))).replace("(", "[").replace(")", "]"))

    series = run2(matters, parts, nums, rank)
    with open(direc+"/Sp4s1a1nf0.txt", 'w+') as f:
        f.write(series)
        f.close()
    print(series)

    end = time.time()
    print(end - start)

    url = os.environ.get("SLACK_WEBHOOK_URL", "")

    title = ("Calculation Complete :zap:")
    message = 'Hilbert series calculation completed\nTime: %d' % (end - start)

    slack_data = {
        "username": "NotificationBot",
        "icon_emoji": ":satellite:",
        "attachments": [
            {
                "color": "#9733EE",
                "fields": [
                    {
                        "title": title,
                        "value": message,
                        "short": "false",
                    }
                ]
            }
        ]
    }

    byte_length = str(sys.getsizeof(slack_data))
    headers = {'Content-Type': "appication.json", 'Content-Length': byte_length}
    response = requests.post(url, data=json.dumps(slack_data), headers=headers)
    if response.status_code != 200:
        raise Exception(response.status_code, response.text)


# if __name__ == '__main__':
#     start = time.time()
#     rank = 6
#     count = [0, 0, 0, 0, 0, 2, 2, 0, 0]
#     counts=str(count).replace("[", "{").replace("]", "}")
#     order = 6
#     nums = num(count)
#
#     matters = matter(rank,count)
#     parts = ast.literal_eval(
#         str(list(product(range(0, order + 1), repeat=len(ast.literal_eval(nums))))).replace("(", "[").replace(")", "]"))
#
#     series = run(matters, parts, nums, rank)
#     with open(direc+"/SU5s2S2nf0_1.txt", 'w+') as f:
#         f.write(series)
#         f.close()
#     print(series)
#
#     end = time.time()
#     print(end - start)
#     text = 'Hilbert series calculation completed\nTime: %d' % (end - start)
#     asyncio.run(tel(text))


# result=ast.literal_eval(res)
# ll=[]
# for i in result:
#     l=len(i)
#     p = ""
#     for j in range(0,l-1):
#         if i[j]!=0:
#             p=p+mark[j]+"^"+str(i[j])
#     ll.append([p,i[l-1]])


# for i in ll:
#     print(i)

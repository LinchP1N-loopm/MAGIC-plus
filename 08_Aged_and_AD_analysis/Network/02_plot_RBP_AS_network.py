from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text

# Parameters
ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'figures'
OUT.mkdir(exist_ok=True)
COLORS = {3:'#D64C59',4:'#008DCB',16:'#8C70A5',8:'#358C84'}
plt.rcParams.update({'font.family':'Arial','font.size':7,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
edges = pd.read_csv(ROOT/'csvs/34_AS_only_module_internal_edges.tsv.gz',sep='\t')
full = pd.read_csv(ROOT/'csvs/30_AS_only_RBP_AS_edges.tsv.gz',sep='\t')
nodes = pd.read_csv(ROOT/'csvs/31_AS_only_module_nodes.tsv',sep='\t')
mp=nodes.set_index('node').module
full['RBP_module']=('RBP|'+full.RBP_name).map(mp)
full['AS_module']=('AS|'+full.EventID).map(mp)
inside=full[full.RBP_module==full.AS_module]
assert set(zip(inside.gene_id,inside.EventID))==set(zip(edges.gene_id,edges.EventID))
full[full.RBP_module!=full.AS_module].to_csv(OUT/'cross_module_edges.tsv',sep='\t',index=False)
# AS_class来自原差异剪接检出分类
class_map={'Shared':'Both comparisons','Both region-specific':'Both comparisons',
           'AD only':'AD vs Young only','Aged only':'Aged vs Young only'}
category=edges.drop_duplicates('EventID').set_index('EventID').AS_class.map(class_map).to_dict()
PALETTE={'AD vs Young only':'#D87830','Aged vs Young only':'#2186B7','Both comparisons':'#8A60A8'}
node_colors={'AS|'+e:PALETTE[c] for e,c in category.items()}
edges['Detection_group']=edges.EventID.map(category)
assert edges.Detection_group.notna().all()
layout_rows=[]

def edge_width(r):
    return .45 + 5*(abs(r)-.5)


def network(ax,module):
    d=edges[edges.module==module]
    g=nx.from_pandas_edgelist(d,'RBP_node','AS_node',['sample_adj_r'])
    event=d[['AS_node','EventID','AS_gene','AS_type']].drop_duplicates('AS_node')
    label=dict(zip(d.RBP_node,d.RBP_name)) | dict(zip(event.AS_node,event.AS_gene+'-'+event.AS_type))
    for node in event.loc[event.AS_gene.eq('Nme3'),'AS_node']: label[node] += '†'
    pos=nx.kamada_kawai_layout(g,weight=None)
    xy=np.array(list(pos.values()));xy-=xy.mean(axis=0)
    _,_,rot=np.linalg.svd(xy,full_matrices=False);xy=xy@rot.T;xy/=np.abs(xy).max(axis=0)
    pos=dict(zip(pos,xy))
    for sign,style in [(1,'solid'),(-1,(0,(3,2)))]:
        ee=[(u,v) for u,v,e in g.edges(data=True) if np.sign(e['sample_adj_r'])==sign]
        nx.draw_networkx_edges(g,pos,edgelist=ee,ax=ax,edge_color='#8C969F',alpha=.85,style=style,
             width=[edge_width(g[u][v]['sample_adj_r']) for u,v in ee])
    nx.draw_networkx_nodes(g,pos,nodelist=d.RBP_node.unique().tolist(),ax=ax,node_size=8,node_color='#41484D',linewidths=0)
    nx.draw_networkx_nodes(g,pos,nodelist=event.AS_node.tolist(),ax=ax,node_size=19,node_color=[node_colors[n] for n in event.AS_node],linewidths=0)
    ax.set_xlim(-1.5,1.5);ax.set_ylim(-1.35,1.35)
    txt=[ax.text(*pos[n],label[n],ha='center',va='center',fontsize=7,fontstyle='italic',color=node_colors[n] if n.startswith('AS|') else '#363B40') for n in g]
    if mask=='recommended' and module==3:
        offsets={'Vdac3-A5SS':(-.22,.35),'Nme3-SE†':(-.25,-.16),'Emc10-SE':(.05,-.25),'Smg7':(-.04,.23),'Zbtb2':(-.16,.21),'Csrp2':(0,.20),'Trim56':(.28,.08),'App-SE':(.29,.02),'Exoc7-SE':(.18,.20),'Icam1':(.15,-.17),'Eif3e':(-.13,-.18),'Mrpl10':(-.02,.23)}
        for n,t in zip(g,txt):
            dx,dy=offsets[label[n]];t.set_position((pos[n][0]+dx,pos[n][1]+dy))
    else:
        np.random.seed(42)
        adjust_text(txt,x=xy[:,0],y=xy[:,1],ax=ax,expand=(1.2,1.5),force_text=(.5,.7),force_static=(.5,.8),
            iter_lim=900,min_arrow_len=4,arrowprops=dict(arrowstyle='-',color='#BDBDBD',lw=.3))
    ax.set_title(f'Module {module}',loc='left',color='#30363B',fontsize=8,pad=7)
    ax.text(1,1.035,f'{d.RBP_node.nunique()} RBPs',transform=ax.transAxes,ha='right',fontsize=6.5,color='#777777')
    ax.grid(False);ax.axis('off')
    for n,t in zip(g,txt):layout_rows.append({'module':module,'node':n,'label':label[n],'x':pos[n][0],'y':pos[n][1],'label_x':t.get_position()[0],'label_y':t.get_position()[1]})


def save(fig,name):
    fig.savefig(OUT/f'{name}.pdf')
    fig.savefig(OUT/f'{name}.svg')
    fig.savefig(OUT/f'{name}.png',dpi=600)
    plt.close(fig)


# 筛选：Top 2 加桥接边；严格 Top 3 供比较
all_edges=edges.copy()
all_edges.to_csv(OUT/'Supplement_all_module_internal_edges.tsv',sep='\t',index=False)
full.to_csv(OUT/'Supplement_all_AS_only_edges.tsv',sep='\t',index=False)
selections=[]
counts=[]
for module in [3,4]:
    d=all_edges[(all_edges.module==module)&(all_edges.sample_adj_padj<.01)].copy()
    d['strength']=d.sample_adj_r.abs()
    d=d.sort_values(['strength','RBP_name','gene_id'],ascending=[False,True,True])
    d['rank_within_event']=d.groupby('EventID').cumcount()+1
    g=nx.from_pandas_edgelist(d,'RBP_node','AS_node',['strength'])
    tree=nx.maximum_spanning_tree(g,weight='strength')
    while True:
        leaves=[n for n in tree if n.startswith('RBP|') and tree.degree(n)<2]
        if not leaves:break
        tree.remove_nodes_from(leaves)
    backbone={frozenset(e) for e in tree.edges}
    d['top2']=d.rank_within_event<=2
    d['bridge_backbone']=[frozenset((r.RBP_node,r.AS_node)) in backbone for r in d.itertuples()]
    d['recommended']=d.top2|d.bridge_backbone
    d['strict_top3']=d.rank_within_event<=3
    selections.append(d)
    for version,mask in [('Full',d.index==d.index),('Top2',d.top2),('Top3',d.strict_top3),('Top2_plus_bridges',d.recommended)]:
        z=d[mask];h=nx.from_pandas_edgelist(z,'RBP_node','AS_node')
        counts.append(dict(module=module,version=version,n_RBP=z.RBP_node.nunique(),n_AS=z.EventID.nunique(),n_edges=len(z),components=nx.number_connected_components(h)))
selected=pd.concat(selections)
selected.to_csv(OUT/'Main_modules_full_edges_with_selection.tsv',sep='\t',index=False)
pd.DataFrame(counts).to_csv(OUT/'Selection_counts.tsv',sep='\t',index=False)
print(pd.DataFrame(counts).to_string(index=False))
handles=[Line2D([],[],marker='o',ls='',color='#41484D',ms=3,label='RBP')]+[
    Line2D([],[],marker='o',ls='',color=c,ms=4,label=k) for k,c in PALETTE.items()]
sign_handles=[Line2D([],[],color='#8C969F',lw=1,label='Positive',ls='-'),
    Line2D([],[],color='#8C969F',lw=1,label='Negative',ls=(0,(3,2)))]+[
    Line2D([],[],color='#8C969F',lw=edge_width(r),label=str(r)) for r in [.55,.65,.75]]
for mask,name,subtitle in [('recommended','Main_Top2_plus_bridges_120mm','Top 2 per AS event + connecting edges; adjusted FDR < 0.01'),('strict_top3','Alternative_strict_Top3_120mm','Top 3 per AS event; adjusted FDR < 0.01')]:
    edges=selected[selected[mask]].copy()
    edges.to_csv(OUT/f'{name}_edges.tsv',sep='\t',index=False)
    fig,axs=plt.subplots(1,2,figsize=(120/25.4,75/25.4))
    fig.subplots_adjust(left=.025,right=.98,bottom=.17,top=.69,wspace=.10)
    for ax,module in zip(axs,[3,4]):network(ax,module)
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,1.005),ncol=2,frameon=False,fontsize=6.5,handlelength=1,columnspacing=1.5,labelspacing=.35)
    fig.legend(handles=sign_handles,loc='upper center',bbox_to_anchor=(.5,.855),ncol=5,frameon=False,fontsize=6,handlelength=1.8,columnspacing=.8)
    fig.text(.5,.767,'Line width: |r_adj|',ha='center',fontsize=6,color='#555555')
    fig.text(.04,.09,subtitle,fontsize=6,color='#555555')
    fig.text(.04,.035,'† Nme3: gene results missing in 2 of 3 AS-significant contexts.',fontsize=6,color='#555555')
    save(fig,name)
pd.DataFrame(layout_rows).to_csv(OUT/'network_layouts.tsv',sep='\t',index=False)

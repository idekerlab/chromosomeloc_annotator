import json
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple
import re
import requests
from ndex2.cx2 import RawCX2NetworkFactory, CX2Network
import ndex2.client


HUMAN_CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


@dataclass
class Config:
    client: Optional[ndex2.client.Ndex2] = None
    chrom_map_path: Optional[str] = None
    interaction_gene_attrs: List[str] = field(default_factory=lambda: ["represents", "name"])
    interaction_uuid_attr: str = "HCX::interactionNetworkUUID"
    interaction_uuid_network_attr: str = "HCX::interactionNetworkUUID"
    gene_list_delim: str = ","
    species: str = "human"
    cache_file: Optional[str] = None
    ndex_server: str = "https://public.ndexbio.org"
    ndex_username: Optional[str] = None
    ndex_password: Optional[str] = None
    allow_missing_uuid: bool = False


def load_json_file(path: str) -> Any:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def load_cx2_from_file(path: str, err_stream=sys.stderr) -> CX2Network:
    """
    Load a CX2 network from a file path (or '-' for stdin) into memory.
    """
    err_stream.write('@@PROGRESS 10\n')
    err_stream.write('@@MESSAGE Loading network\n')
    err_stream.flush()

    if path == "-":
        data = json.load(sys.stdin)
    else:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
    cx2_netfactory = RawCX2NetworkFactory()
    return cx2_netfactory.get_cx2network(data)


def _normalize_gene(gene: str) -> str:
    return gene.split(".")[0].strip().upper()


def get_chromosome_map() -> Dict[str, str]:
    """
    Parses the non_alt_loci_set.json file from HGNC to build a gene symbol → chromosome map for non-alt loci genes. This is used as a fallback when no local map file is provided, and is also the source of truth for which genes are considered non-alt loci for the purposes of annotation.
    """
    data = load_json_file(str(Path(__file__).with_name("non_alt_loci_set.json")))
    for entry in data['response']['docs']:
        gene = _normalize_gene(entry['symbol'])
        if 'location' not in entry:
            chrom= "Un"
        else:
            chrom = re.sub('[p|q].*$', '', entry['location'].split(' ')[0])
            if chrom.startswith("mito"):
                chrom = "M"
            elif chrom.startswith("reserved") or chrom.startswith("not"):
                chrom = "Un"
        
        yield gene, f"chr{chrom}"


def _apply_pie_style(data: Any, chromosomes: List[str]) -> None:
    pie_attr_list = ",".join([f"{c}_count" for c in chromosomes])
    colors = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
        "#7f7f7f", "#bcbd22", "#17becf", "#393b79", "#637939", "#8c6d31", "#843c39",
        "#7b4173", "#3182bd", "#e6550d", "#31a354", "#756bb1", "#636363", "#969696",
        "#bdbdbd", "#d9d9d9", "#9e9ac8", "#bcbddc"
    ]
    color_list = ",".join(colors[: len(chromosomes)])
    graphic = f'org.cytoscape.PieChart: attributelist="{pie_attr_list}" colorlist="{color_list}"'

    vp = get_aspect(data, "visualProperties")
    if vp is None:
        vp = [{"default": {"node": {"NODE_CUSTOMGRAPHICS_1": graphic}}, "nodeMapping": {}, "edgeMapping": {}, "network": {}}]
        if isinstance(data, list):
            data.append({"visualProperties": vp})
        else:
            data["visualProperties"] = vp
        return
    # Add to first style default node properties without clobbering other fields
    first_style = vp[0]
    defaults = first_style.setdefault("default", {})
    node_defaults = defaults.setdefault("node", {})
    node_defaults["NODE_CUSTOMGRAPHICS_1"] = graphic


def get_network_from_ndex(uuid: str, config: Config) -> CX2Network:
    factory = RawCX2NetworkFactory()
    return factory.get_cx2network(config.client.get_network_as_cx2_stream(uuid).json())


def get_hierarchy_interactome(net: CX2Network, config: Config) -> CX2Network:
    """
    Get the hierarchy interactome for the given network and configuration.
    using ndex_client to query for the hierarchy based on uuid and return
    the CX2Network object
    """
    if config.interaction_uuid_network_attr not in net.get_network_attributes():
        None
    hierarchy_uuid = net.get_network_attributes()[config.interaction_uuid_network_attr]
    return get_network_from_ndex(hierarchy_uuid, config)


def get_node_interactome(node: Dict[str, Any], config: Config) -> Optional[CX2Network]:
    """
    Get the interactome for a given node based on the interaction UUID attribute and configuration. Returns None if the attribute is missing or if allow_missing_uuid is True and the attribute is missing.
    """
    if config.interaction_uuid_attr not in node.get("attributes", {}):
        return None
        
    interaction_uuid = node["attributes"][config.interaction_uuid_attr]
    return get_network_from_ndex(interaction_uuid, config)

def get_node_id_to_name_map(net: CX2Network) -> Dict[str, str]:
    """
    Build a mapping of node ID to node name for a given network. 
    This is used to link hierarchy nodes to interactome nodes based on shared names.
    """
    name_id_map = {}
    for node_id, node_obj in net.get_nodes().items():
        if "name" in node_obj['v']:
            name_id_map[node_id] = node_obj["v"]["name"]
    return name_id_map


def annotate_network(
    net: CX2Network,
    config: Config,
    err_stream=sys.stderr,
) -> Any:
    err_stream.write("@@PROGRESS 30\n")
    err_stream.write("@@MESSAGE Building chromosome map\n")
    chrom_map = dict(get_chromosome_map())

    chromosomes = sorted(set(chrom_map.values()))

    err_stream.write("@PROGRESS 50\n")
    err_stream.write("@@MESSAGE Annotating nodes\n")
    hierarchy_net = get_hierarchy_interactome(net, config)
    hier_node_map = get_node_id_to_name_map(hierarchy_net)
    for node_id, node_obj in net.get_nodes().items():
        chrom_count = Counter()
        if 'HCX::members' in node_obj['v']:
            members = node_obj['v']['HCX::members']
        elif 'HCX::memberNames' in node_obj['v']:
            members = node_obj['v']['HCX::memberNames'].split(',')
        else:
            continue
        interactome_net = get_node_interactome(node_obj, config)
        if interactome_net is not None:
            the_node_map = get_node_id_to_name_map(interactome_net)
        else:
            the_node_map = hier_node_map
        for member_id in members:
            if not isinstance(member_id, str):
                member_name = member_id
            else:
                member_name = the_node_map.get(member_id)
            if member_name is None:
                chrom_count["chrUn"] += 1
                continue
            if member_name in chrom_map:
                chrom = chrom_map[member_name]
                chrom_count[chrom] += 1
            else:
                chrom_count["chrUn"] += 1
        # update chrom_counts for this node        
            for chrom in chromosomes:
                if chrom in chrom_count:
                    net.add_node_attribute(node_id, f"{chrom}",
                                           chrom_count[chrom], datatype='integer')
                else:
                    net.add_node_attribute(node_id, f"{chrom}", 0, datatype='integer')

        # err_stream.write('node: {}\n'.format(node_id) + ' chrom_count: {}\n'.format(chrom_count))

                        
    err_stream.write("@@PROGRESS 90\n")
    err_stream.write("@@MESSAGE Formatting updated network\n")
    return [{
        "action": "addNetworks",
        "data": [net.to_cx2()]
    }]


def run_update(
    input_source: str,
    output_stream,
    config: Config,
    updatedby: str = "chromloc",
) -> None:
    net = load_cx2_from_file(input_source)
    result = annotate_network(net, config)
    json.dump(result, output_stream, indent=2)
    output_stream.write("\n")
    output_stream.flush()

import re
def _check_list_all_in(l, config):
    for k in l:
        if k not in config:
            return False
    return True
def _decide_mode(l0, l1, config):
    if _check_list_all_in(l0, config) is True:
        return 0
    else:
        if _check_list_all_in(l1, config) is True:
            return 1
        else:
            raise ValueError('Missing columns to identify effect size and SE from GWAS. We need either {l0} or {l1}.')
def _fill_config_w_default_if_needed(config, default_dict):
    entries = list(default_dict.keys())
    for entry in entries:
        if entry not in config:
            config[entry] = default_dict[entry]
def _try_to_format(pattern, key_value_pairs):
    for k, v in key_value_pairs.items():
        kpattern = '{' + k + '}'
        if kpattern in pattern:
            pattern = re.sub(kpattern, v, pattern)
    return pattern
def _try_fill_config(config, key, value, value_pools=None):
    if key not in config:
        config[key] = value
    if value_pools is not None:
        _check_val_in_pool(config[key], value_pools)
def _check_val_in_pool(val, pool):
    if val not in pool:
        raise ValueError(f'{val} not in {pool}.')
def _check_desired_wildcards(str_, wc_list):
    for wc in wc_list:
        if wc not in str_:
            raise ValueError(f'Expect {wc} but not found.')
def _try_format_for_list(pat, dict_, l, key):
    o = []
    for i in l:
        dict_[key] = str(i)
        o.append(_try_to_format(pat, dict_))
    return o
import re, random, itertools

def sorted_nicely(l):
    """ Sorts the given iterable in the way that is expected.
 
    Required arguments:
    l -- The iterable to be sorted.
 
    """
    convert = lambda text: int(text) if text.isdigit() else text
    #Remove the .lower() if you want uppercase letters to come before all lowercase letters
    alphanum_key = lambda key: [convert(c.lower()) for c in re.split('([0-9]+)', key)]
    return sorted(l, key = alphanum_key)

def concat_str_list(str_list, delimiter=' '):
    concat_str = ''
    for elem in str_list:
        concat_str += elem + delimiter
    return concat_str

def indices(lst, val):
    '''
    lst: list
        list to search
    val: *
        value to search the list for

    return: list
        All indices matching val. If there are no indices
        in lst matching val, return empty list
    Purpose: output indices where the val is found in lst. If val is
        not found in lst, output empty list
    '''
    return [i for i, item in enumerate(lst) if item is val or item == val]

def random_index(lst):
    '''
    lst: list
        List
    
    return: int
        random number that could be an index

    Purpose: get a valid random index
    '''
    if len(lst) == 0:
        raise ValueError('lst has len 0.')
    return random.randint(0, len(lst) - 1)

def random_value(lst):
    '''
    lst: list
        List
    
    return: int
        random element of lst

    Purpose: return one of the values in lst at random
    '''
    if len(lst) == 0:
        raise ValueError('lst has len 0.')
    return lst[random_index(lst)]

def flatten_list(lst):
    '''
    lst: list
        List of lists

    return: list
        1D list which is the flattened list
    Purpose: Flatten a list of lists
    '''
    return list(itertools.chain(*lst))

def sort_list_by_col(lst, col, reverse=False):
    return sorted(lst,key=lambda l:l[col], reverse=reverse)


def multi_delete(list_, indices_to_delete):
    '''
    list_: iterable
        list or iterable
    indices_to_delete: iterable
        indices of list_ to delete from list_
    
    return: list
       list of entries in iterable but now with indices_to_delete deleted
    Purpose: delete many indices from an iterable at once to avoid 
        clashes of indices changing as elements are deleted
    '''
    indexes = sorted(list(indices_to_delete), reverse=True)
    for index in indexes:
        del list_[index]
    return list_


def randomsubset_not_in_other_set(n_choose, subset_to_avoid, total_set):
    '''
    n_choose: int
        number of indices to choose

    subset_to_avoid: iterable
        set of indices that you don't want to include in your subset

    total_set: iterable or None
        The set of indices to pick from. If None, then

    Return: np array, shape: ,n_choose
        A list of indices to include in the subset.

    Purpose: Determine a random subset of the entire population of indices
        (total_set) to include
       in a subset.
    '''
    if type(subset_to_avoid) is not set:
        subset_to_avoid = set(subset_to_avoid)
    if type(total_set) is not set:
        total_set = set(total_set)

    valid_set_to_choose_from = total_set.difference(subset_to_avoid)
    return np.array(random.sample(valid_set_to_choose_from, n_choose))


def sort_by_col(data, col):
    '''
    data: numpy array, shape (at least one column)
        array to sort

    col: int
        column index to sort

    return: the data sorted by the specified column

    purpose: Sort data by the specified column
    '''
    try:
        sorted_data = data[np.argsort(data[:,col])]
    except:
        return None
    return sorted_data


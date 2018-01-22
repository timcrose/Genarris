def check_int(arg, arg_name):
    if not isinstance(arg, int):
        raise TypeError(arg_name + " must be an integer.")

def check_int_or_none(arg, arg_name):
    if not isinstance(arg, int) and not arg is None:
        raise TypeError(arg_name + " must be an integer or None.")

def check_float(arg, arg_name):
    if not isinstance(arg, (float, int)):
        raise TypeError(arg_name + " must be a float or integer.")

def check_float_or_none(arg, arg_name):
    if not isinstance(arg, (float, int)) and not arg is None:
        raise TypeError(arg_name + " must be a float or integer.")

def check_list(arg, arg_name):
    if not isinstance(arg, list):
        raise TypeError(arg_name + " must be a list.")

def check_list_or_none(arg, arg_name):
    if not isinstance(arg, list) and not arg is None:
        raise TypeError(arg_name + " must be a list or None.")

def check_dict(arg, arg_name):
    if not isinstance(arg, dict):
        raise TypeError(arg_name + " must be a dict.")

def check_dict_or_none(arg, arg_name):
    if not isinstance(arg, dict) and not arg is None:
        raise TypeError(arg_name + " must be a dict or None.")

def check_type(arg, arg_name, type_list):
    if not isinstance(arg, type_list):
        raise TypeError(arg_name + " must be of one of: " + str(type_list)) 

def check_type_or_none(arg, arg_name, type_list):
    if not isinstance(arg, type_list) and not arg is None:
        raise TypeError(arg_name + " must be of one of: " + str(type_list) +
                " or None")

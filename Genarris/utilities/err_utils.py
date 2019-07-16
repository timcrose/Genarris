
def handle_error(e=None, err_message='Error', fail_gracefully=False, 
                verbose=False):
                
    if fail_gracefully:
        if verbose:
            print(err_message)
    else:
        if e is None:
            raise Exception(err_message)
        else:
            if verbose:
                print(err_message)
            raise e

def check_input_vars(vars_type_list, fail_gracefully=False, verbose=False):
    if not (hasattr(vars_type_list, '__iter__') and \
            type(vars_list) is not str and \
            type(vars_list) is not set and \
            type(vars_list) is not dict):
        err_message='vars_type_list must be iterable and not str, set, or ' + \
                    'dict because set and dict are not guaranteed to be in ' + \
                    'the appropriate order. Got vars_type_list: ' + \
                    str(type(vars_type_list))

        handle_error(e=None, err_message=err_message, fail_gracefully=False,
                    verbose=False)

        return None


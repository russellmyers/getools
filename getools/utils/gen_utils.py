class SerialiserMixin:

    def _to_attr_dict(self):

        attr_dict = vars(self)
        attr_dict_out = attr_dict.copy()
        for attr in attr_dict:
            if isinstance(attr_dict[attr],list):
                list_out = []
                for item in attr_dict[attr]:
                    if isinstance(item,SerialiserMixin):
                        list_out.append(item._to_attr_dict())
                    else:
                        list_out.append(item)
                attr_dict_out[attr] = list_out
            elif isinstance(attr_dict[attr], tuple):
                list_out = []
                for item in attr_dict[attr]:
                    if isinstance(item,SerialiserMixin):
                        list_out.append(item._to_attr_dict())
                    else:
                        list_out.append(item)
                attr_dict_out[attr] = tuple(list_out)
            elif isinstance(attr_dict[attr],SerialiserMixin):
                attr_dict_out[attr] = attr_dict[attr]._to_attr_dict()
            else:
                attr_dict_out[attr] = attr_dict[attr]

        return attr_dict_out


    @staticmethod
    def _from_attr_dict(attr_dict):
        obj = SerialiserMixin()
        return obj

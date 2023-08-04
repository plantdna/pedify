
def rgb_to_hex(rgb:list)->str:
    color = '#'
    for i in rgb:
        num = int(i)
        color += str(hex(num))[-2:].replace('x', '0').upper()
    return color


def hex_to_rgb(hex:str)->list:
    r = int(hex[1:3],16)
    g = int(hex[3:5],16)
    b = int(hex[5:7], 16)
    rgb = [r, g, b]
    return rgb
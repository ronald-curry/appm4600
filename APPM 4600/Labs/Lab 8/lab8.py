def linear_interp (x1,x2,y1,y2,alpha):
    y_alpha=y1+(alpha-x1)*(y2-y1)/(x2-x1)
    return y_alpha
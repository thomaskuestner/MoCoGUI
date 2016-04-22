function imview(varargin)
if isempty(varargin)
    nextview
else
    N=length(varargin);
    if N==1&iscell(varargin{1})
        varargin=varargin{1};
        N=length(varargin);
    end
    clf
    h=image(varargin{1});
    colormap(gray(256))
    axis image,axis off
    if N>=2
        set(h,'userdata',{0 1 varargin 0},'ButtonDownFcn','imview')
    end
end

function nextview
h=findobj(gca,'type','image');
a=get(h,'userdata');
N=length(a{3});
if length(a{3}{N})==1
    K=round(a{3}{N});
    N=N-1;
end

point1=round(get(gca,'CurrentPoint'));
point1=[point1(1,2) point1(1,1)];
tic,rbbox;dt=toc;
point2=round(get(gca,'CurrentPoint'));
point2=[point2(1,2) point2(1,1)];

if dt<0.2
    a=get(h,'userdata');
    if a{4}
        a{4}=0;
        set(h,'userdata',a)
        encore=a{4};
    else
        a{4}=1;
        set(h,'userdata',a)
        encore=a{4};
    end
    while encore
        try
            a=get(h,'userdata');
            a{2}=(1-2*(a{1}+a{2}==0|a{1}+a{2}==N+1))*a{2};
            a{1}=a{1}+a{2};
            set(h,'cdata',a{3}{a{1}},'userdata',a)
            encore=a{4};
            pause(0.1)
        catch
            encore=0;
        end
    end
else
    if point1==point2
        if ~exist('K','var')
            K=2;
        end
        point1=max(point1-K,1);
        point2=min(point2+K,size(a{3}{1}));
    end

    for k=1:N
        a{3}{k}=a{3}{k}(point1(1):point2(1),point1(2):point2(2));
    end
    figure
    h=image(a{3}{1});
    colormap(gray(256))
    axis image,axis off
    if N>=2
        set(h,'userdata',a,'ButtonDownFcn','imview')
    end

    if a{4}
        a{4}=0;
        set(h,'userdata',a)
        encore=a{4};
    else
        a{4}=1;
        set(h,'userdata',a)
        encore=a{4};
    end
    while encore
        try
            a=get(h,'userdata');
            a{2}=(1-2*(a{1}+a{2}==0|a{1}+a{2}==N+1))*a{2};
            a{1}=a{1}+a{2};
            set(h,'cdata',a{3}{a{1}},'userdata',a)
            encore=a{4};
            pause(0.1)
        catch
            encore=0;
        end
    end
end

function nextview2
h=findobj(gca,'type','image');
a=get(h,'userdata');
if a{4}
    a{4}=0;
    set(h,'userdata',a)
    encore=a{4};
else
    a{4}=1;
    set(h,'userdata',a)
    encore=a{4};
end
while encore
    try
        a=get(h,'userdata');
        N=length(a{3});
        a{2}=(1-2*(a{1}+a{2}==0|a{1}+a{2}==N+1))*a{2};
        a{1}=a{1}+a{2};
        set(h,'cdata',a{3}{a{1}},'userdata',a)
        encore=a{4};
        pause(0.1)
    catch
        encore=0;
    end
end

function obj3=join_obj(obj1,obj2)
    fldn1  = fieldnames(obj1);
    fldn2  = fieldnames(obj2);
    for i=1:length(fldn1)
        if ~isfield(obj2,fldn1{i})
            emp=repmat([obj1.(fldn1{i})],0); %hack to make an empty entry in field that is not mutual to the objects. Crucially, empty entry will have the correct data type.
            [obj2(:).(fldn1{i})]=deal(emp);
        end
    end
    for i=1:length(fldn2)
        if ~isfield(obj1,fldn2{i})
            emp=repmat([obj2.(fldn2{i})],0);
            [obj1(:).(fldn2{i})]=deal(emp);
        end
    end
    obj3=[obj1,obj2];
end
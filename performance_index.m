function [accuracy, RI, NMI,FMeasure]=performance_index(real_label,our_id)

[accuracy, label_new]=label_map(real_label,our_id);
[~,RI]=RandIndex(real_label,label_new);
NMI=nmi(real_label',label_new');
[precision,recall,F,FMeasure,Accuracy] = calaccuracy(real_label,label_new);

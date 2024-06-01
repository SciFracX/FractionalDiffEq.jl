function Base.summary(io::IO, prob::AbstractFDEProblem)
    type_color, no_color = SciMLBase.get_colorizers(io)
    print(io,
        type_color,
        nameof(typeof(prob)),
        no_color,
        " with uType ",
        type_color,
        typeof(prob.u0),
        no_color,
        " and tType ",
        type_color,
        prob.tspan isa Function ? "Unknown" :
        (prob.tspan === nothing ? "Nothing" : typeof(prob.tspan[1])),
        no_color,
        ". In-place: ",
        type_color,
        SciMLBase.isinplace(prob),
        no_color)
end

function Base.show(io::IO, mime::MIME"text\plain", A::FDDEProblem)
    summary(io, A)
    println(io)
    print(io, "timespan: ")
    show(io, mime, A.tspan)
    println(io)
    print(io, "order: ")
    show(io, mime, A.order)
    println(io)
    print(io, "u0: ")
    show(io, mime, A.u0)
end

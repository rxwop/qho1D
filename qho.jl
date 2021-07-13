using Plots, QuadGK, ColorSchemes, ForwardDiff

# Vector e.g. [1, 4, 5] or collected Range e.g. collect(1:5)
inputstates = [0, 2, 3]
sort!(inputstates)

function constant(n)
    1/sqrt(size(inputstates, 1))
    #2/(n*pi)*(1-cos(n*pi/2))
end

function calculate(states)
    dx = 0.25
    dt = 1/30
    elapse = 10
    m = 0.2
    w = 1

    h = 1.054571817
    A = sqrt(2h*(states[end]+0.5)/w/m)+4
    xran = (-A):dx:A
    yran = 0:dx:elapse
    clrs = ColorSchemes.GnBu_9
    grad1 = cgrad(ColorSchemes.rainbow)
    

    Cₙ = [constant(i) for i in states]

    function GenerateHermite(n)
        Hermite=Function[]
    
        push!(Hermite,x->1);
        push!(Hermite,x->2*x);
    
        for ni in 3:n
            push!(Hermite,x->2*x*Hermite[ni-1](x)-2*(ni-1)*Hermite[ni-2](x))
        end
        return Hermite
    end

    Hermite=GenerateHermite(states[end]+1)

    Ψₙ(x) = [(m*w/pi/h)^(1/4)*Hermite[i+1](sqrt(m*w/h)*x)*exp(-m*w*x^2/2/h)/sqrt(2^i*factorial(big(i))) for i in states]
    dΨₙ(x) = ForwardDiff.derivative(Ψₙ, x)

    φₙ(t) = [exp(-im*w*t*(i+0.5)) for i in states]

    radii = [[0,abs(i)] for i in Cₙ]

    E(n)=h*w*(n+0.5)

    Eₙ = [E(i) for i in states]

    

    anim = @animate for ti in 0:dt:elapse

            Ψ(x) = sum(Ψₙ(x).*φₙ(ti).*Cₙ)
            dΨ(x) = sum(dΨₙ(x).*φₙ(ti).*Cₙ)
            
            ref = real.(Ψ.(xran))
            imf = imag.(Ψ.(xran))
            arg = vcat(atan.(ref,imf), pi, -pi)
            ρ(x) = abs2(Ψ(x))
            V(x) = m*w^2*x^2/2
            J(x) = real(im*h/(2*m)*(Ψ(x)*conj(dΨ(x))-conj(Ψ(x))*dΨ(x)))
            xΨ(x) = x*ρ(x)
            x²Ψ(x) = x^2*ρ(x)
            #expcX = first(quadgk(xΨ, -Inf, Inf))
            #expcX² = first(quadgk(x²Ψ, -Inf, Inf))
            #σₓ = sqrt(expcX²-expcX^2)

            theta = [[0,-i] for i in angle.(φₙ(ti))]

            p1 = plot(xran, imf, ref, ylims = [-1,1], zlims = [-1,1], lc = grad1, labels = false, title = "Ψ(x, t)", zlabel = "Re(Ψ)", ylabel = "Im(Ψ)", line_z = arg, legend = false, camera = [20sin(ti)+55, 25])
            p2 = plot(xran, ρ, ylims = [0,1], color = clrs[6], labels = false, title = "∣Ψ(x, t)∣²", fill = true)
            p3 = plot(theta, radii, proj =:polar, lims = [0,1], showaxis = false, labels = false, palette = [clrs[4], clrs[5], clrs[6], clrs[7], clrs[8]], title = "φ(t)")
            #scatter!(p2, [expcX], [0.5], color = clrs[5], xerror = σₓ, markerstrokecolor=clrs[8], labels = false)
            plot!(p1, xran, zeros(size(xran,1),1), zeros(size(xran,1),1), label = false, alpha = 0.1, color = :black)
            plot!(p2, V, label = "V(x)", color = clrs[5], linestyle=:dash)

            p5 = plot(xran, J, ylims = [-A, A], labels = false, title = "J(x,t)", color = clrs[6])

            l1 = @layout[a b; c{0.4h} d]

            plot(p1, p2, p3, p5, fontfamily = :serif, layout = l1)

    end

    Eplot = vcat(E(states[begin]-2), E(states[begin]-1), Eₙ, E(states[end]+1), E(states[end]+2))
    Cplot = vcat(0, 0, abs2.(Cₙ), 0, 0)
    Psi(x, t) = abs2(sum(Ψₙ(x).*φₙ(t).*Cₙ))
    l2 = @layout[a{0.3w} [b; c]]

    p4 = plot(Cplot, Eplot, xlims = [0, 1], markershapes = :circle, ylabel = "E", xlabel = "∣C∣²", labels = false, line = clrs[4], markercolors = clrs[4])
    hline!(p4, [sum(abs2.(Cₙ).*Eₙ)], color = clrs[9], label = false)
    p8 = surface(xran, yran, Psi, camera = [15, 75], seriescolor = cgrad(clrs), ylabel = "t", zlabel = "∣Ψ(x, t)∣²")
    p6 = contourf(xran, yran, Psi, seriescolor = cgrad(clrs), ylabel = "t")

    stats = plot(p4, p8, p6, layout = l2, fontfamily = :serif)

    savefig(stats, "qho1_$inputstates _stats.png")

    gif(anim, "qho1_$inputstates.gif", fps = 1/dt)
end

@time calculate(inputstates)

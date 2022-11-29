(pwd() != @__DIR__) && cd(@__DIR__) # allow starting app from bin/ dir

module EmbedViewer

using Genie

const up = Genie.up
export up

function main()
  Genie.genie(; context = @__MODULE__)
end

end

const UserApp = EmbedViewer
EmbedViewer.main()
